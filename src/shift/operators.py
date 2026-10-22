"""Auditable numerical differential operators for structured Cartesian data."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


def _coordinate(value: ArrayLike, name: str, *, minimum_points: int) -> FloatArray:
    coordinate = np.asarray(value, dtype=float)
    if coordinate.ndim != 1 or coordinate.size < minimum_points:
        raise ValueError(
            f"{name} must be a 1-D coordinate with at least {minimum_points} points"
        )
    if not np.all(np.isfinite(coordinate)) or not np.all(np.diff(coordinate) > 0.0):
        raise ValueError(f"{name} must be finite and strictly increasing")
    return coordinate


def _finite_difference_weights(
    points: FloatArray, evaluation_point: float, derivative_order: int
) -> FloatArray:
    """Polynomial finite-difference weights on one local nonuniform stencil."""
    if derivative_order < 0 or derivative_order >= points.size:
        raise ValueError("derivative order must be smaller than stencil size")
    offsets = points - evaluation_point
    scale = float(np.max(np.abs(offsets)))
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("finite-difference stencil must contain distinct points")
    scaled = offsets / scale
    matrix = np.vstack(
        [scaled**power / math.factorial(power) for power in range(points.size)]
    )
    rhs = np.zeros(points.size, dtype=float)
    rhs[derivative_order] = 1.0
    weights = np.linalg.solve(matrix, rhs)
    return weights / scale**derivative_order


def _stencil_slice(index: int, point_count: int, stencil_size: int) -> slice:
    left = stencil_size // 2
    start = index - left
    start = max(0, min(start, point_count - stencil_size))
    return slice(start, start + stencil_size)


def differentiate(
    values: ArrayLike,
    coordinate: ArrayLike,
    *,
    axis: int,
    derivative_order: int,
) -> FloatArray:
    """Differentiate on a nonuniform 1-D axis using local polynomial stencils.

    First derivatives use three points and are second-order accurate on smooth
    systematically refined grids. Second derivatives use four points and are
    second-order accurate, including one-sided boundary stencils.
    """
    if derivative_order not in (1, 2):
        raise ValueError("SHIFT Phase 2 supports derivative_order 1 or 2")
    stencil_size = 3 if derivative_order == 1 else 4
    x = _coordinate(coordinate, "coordinate", minimum_points=stencil_size)
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError("values must contain only finite samples")
    normalized_axis = axis if axis >= 0 else array.ndim + axis
    if normalized_axis < 0 or normalized_axis >= array.ndim:
        raise ValueError("axis is out of bounds for values")
    if array.shape[normalized_axis] != x.size:
        raise ValueError("coordinate length must match the differentiated axis")

    moved = np.moveaxis(array, normalized_axis, 0)
    result = np.empty_like(moved, dtype=float)
    for index in range(x.size):
        stencil = _stencil_slice(index, x.size, stencil_size)
        points = x[stencil]
        weights = _finite_difference_weights(
            points, float(x[index]), derivative_order
        )
        result[index] = np.tensordot(weights, moved[stencil], axes=(0, 0))
    return np.moveaxis(result, 0, normalized_axis)


def first_derivative(
    values: ArrayLike,
    coordinate: ArrayLike,
    *,
    axis: int,
) -> FloatArray:
    """Second-order first derivative on a finite nonuniform coordinate."""
    return differentiate(values, coordinate, axis=axis, derivative_order=1)


def second_derivative(
    values: ArrayLike,
    coordinate: ArrayLike,
    *,
    axis: int,
) -> FloatArray:
    """Second-order second derivative on a finite nonuniform coordinate."""
    return differentiate(values, coordinate, axis=axis, derivative_order=2)


@dataclass(frozen=True, slots=True)
class CartesianGrid:
    """Explicit structured Cartesian coordinates in metres."""

    x_m: FloatArray
    y_m: FloatArray
    z_m: FloatArray

    def __post_init__(self) -> None:
        object.__setattr__(self, "x_m", _coordinate(self.x_m, "x_m", minimum_points=4))
        object.__setattr__(self, "y_m", _coordinate(self.y_m, "y_m", minimum_points=4))
        object.__setattr__(self, "z_m", _coordinate(self.z_m, "z_m", minimum_points=4))

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.x_m.size, self.y_m.size, self.z_m.size

    @property
    def maximum_spacing_m(self) -> float:
        return float(
            max(
                np.max(np.diff(self.x_m)),
                np.max(np.diff(self.y_m)),
                np.max(np.diff(self.z_m)),
            )
        )


def _validate_scalar_support(field: ArrayLike, grid: CartesianGrid) -> FloatArray:
    array = np.asarray(field, dtype=float)
    if array.ndim < 3 or array.shape[-3:] != grid.shape:
        raise ValueError("scalar field must have trailing spatial shape (nx, ny, nz)")
    if not np.all(np.isfinite(array)):
        raise ValueError("scalar field must contain only finite values")
    return array


def _validate_vector_support(field: ArrayLike, grid: CartesianGrid) -> FloatArray:
    array = np.asarray(field, dtype=float)
    if (
        array.ndim < 4
        or array.shape[-4:-1] != grid.shape
        or array.shape[-1] != 3
    ):
        raise ValueError("vector field must have trailing shape (nx, ny, nz, 3)")
    if not np.all(np.isfinite(array)):
        raise ValueError("vector field must contain only finite values")
    return array


def scalar_gradient(field: ArrayLike, grid: CartesianGrid) -> FloatArray:
    """Gradient of a scalar field; output has trailing vector dimension 3."""
    value = _validate_scalar_support(field, grid)
    axes = (-3, -2, -1)
    coordinates = (grid.x_m, grid.y_m, grid.z_m)
    derivatives = [
        first_derivative(value, coordinate, axis=axis)
        for axis, coordinate in zip(axes, coordinates, strict=True)
    ]
    return np.stack(derivatives, axis=-1)


def scalar_laplacian(field: ArrayLike, grid: CartesianGrid) -> FloatArray:
    """Cartesian scalar Laplacian on the structured grid."""
    value = _validate_scalar_support(field, grid)
    axes = (-3, -2, -1)
    coordinates = (grid.x_m, grid.y_m, grid.z_m)
    result = np.zeros_like(value, dtype=float)
    for axis, coordinate in zip(axes, coordinates, strict=True):
        result += second_derivative(value, coordinate, axis=axis)
    return result


def vector_gradient(field: ArrayLike, grid: CartesianGrid) -> FloatArray:
    """Velocity gradient A_ij = partial u_i / partial x_j."""
    value = _validate_vector_support(field, grid)
    component_gradients = [
        scalar_gradient(value[..., component], grid) for component in range(3)
    ]
    return np.stack(component_gradients, axis=-2)


def vector_laplacian(field: ArrayLike, grid: CartesianGrid) -> FloatArray:
    """Componentwise Cartesian vector Laplacian."""
    value = _validate_vector_support(field, grid)
    return np.stack(
        [scalar_laplacian(value[..., component], grid) for component in range(3)],
        axis=-1,
    )


def divergence(field: ArrayLike, grid: CartesianGrid) -> FloatArray:
    """Divergence of a vector field."""
    gradient = vector_gradient(field, grid)
    return np.trace(gradient, axis1=-2, axis2=-1)


def curl(field: ArrayLike, grid: CartesianGrid) -> FloatArray:
    """Curl of a vector field."""
    gradient = vector_gradient(field, grid)
    return np.stack(
        (
            gradient[..., 2, 1] - gradient[..., 1, 2],
            gradient[..., 0, 2] - gradient[..., 2, 0],
            gradient[..., 1, 0] - gradient[..., 0, 1],
        ),
        axis=-1,
    )


def temporal_derivative(
    values: ArrayLike,
    time_s: ArrayLike,
    *,
    axis: int = 0,
) -> FloatArray:
    """Second-order temporal derivative on a nonuniform time coordinate."""
    return first_derivative(values, time_s, axis=axis)


def wall_normal_derivative(
    samples: ArrayLike,
    wall_distance_m: ArrayLike,
    *,
    axis: int,
) -> FloatArray:
    """One-sided derivative at the wall from explicit wall-normal samples.

    The first coordinate is the wall and must be zero. Distances increase into
    the fluid. The returned value is the derivative at that first sample.
    """
    distance = _coordinate(wall_distance_m, "wall_distance_m", minimum_points=3)
    if not np.isclose(distance[0], 0.0, rtol=0.0, atol=1e-15):
        raise ValueError("wall_distance_m must start at zero at the wall")
    derivative = first_derivative(samples, distance, axis=axis)
    normalized_axis = axis if axis >= 0 else derivative.ndim + axis
    if normalized_axis < 0 or normalized_axis >= derivative.ndim:
        raise ValueError("axis is out of bounds for differentiated result")
    return np.take(derivative, 0, axis=normalized_axis)
