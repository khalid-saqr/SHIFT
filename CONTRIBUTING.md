# Contributing to SHIFT

SHIFT is scientific software. Correctness, provenance, and interpretability take priority over feature count.

## Before contributing

Read:

1. `docs/SCIENTIFIC_CONTRACT.md`
2. `docs/EXPERIMENT_PROTOCOL.md`
3. `docs/ARCHITECTURE.md`

The scientific contract is normative for the source04 implementation.

## Change classes

### Scientific-contract changes

A change to a definition, sign convention, state coordinate, distance, admissible flow class, or closure proposition requires:

- an explicit issue describing the theoretical reason;
- identification of the affected theory version;
- tests showing the old and new behavior;
- a documented migration path;
- no motivation based solely on obtaining a desired experimental outcome.

### Numerical-method changes

A new derivative, interpolation, filtering, solver, pressure reconstruction, wall-gradient, or mapping method requires:

- a stated mathematical operator;
- unit tests on fields with known solutions;
- convergence or repeatability evidence appropriate to the operator;
- preservation of dimensional units until explicit normalization;
- uncertainty reporting that can feed threshold derivation.

### Data-adapter changes

An adapter for CFD, PIV/PTV, 4D-flow MRI, or another source must declare what is directly observed, reconstructed, projected, filtered, or unavailable. Missing state variables must not be silently manufactured.

## Documentation math

All Markdown documentation uses GitHub MathJax with dollar delimiters only:

- inline mathematics: `$...$`;
- displayed mathematics: `$$...$$`, with the delimiters on their own lines for multiline equations.

Backslash-style TeX math delimiters are not permitted in Markdown. CI enforces this convention across every `*.md` file in the repository.

## Development workflow

```bash
python -m pip install -e '.[dev]'
pytest
ruff check .
```

Every pull request should be small enough to review scientifically as well as programmatically.

## Pull-request checklist

- [ ] Scientific definitions are unchanged, or the theory-version change is explicit.
- [ ] New public behavior has tests.
- [ ] Units and coordinate conventions are documented.
- [ ] Markdown equations use dollar-delimited MathJax.
- [ ] Numerical thresholds have stated provenance.
- [ ] Undefined quantities remain undefined rather than being replaced by zero.
- [ ] No closure outcome was used to tune a tolerance, weight, mapping, or pair-selection rule.
- [ ] Documentation matches the implementation.
