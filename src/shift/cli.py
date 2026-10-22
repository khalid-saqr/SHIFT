"""Public command-line entry point for SHIFT."""

import argparse
import json
import platform
import sys

from . import __version__

PROJECT_VERSION = __version__
CONTRACT_SUMMARY = (
    "Full vector WSS -> independently qualified near-wall source04 state; "
    "thresholds are frozen from uncertainty before closure outcomes."
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="shift",
        description="Shear-equivalent Hemodynamic Invariance Falsification Test",
    )
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("info", help="show project identity and scientific contract summary")
    sub.add_parser("doctor", help="show the local runtime needed for bug reports")

    benchmark = sub.add_parser(
        "benchmark", help="run a locked analytical regression benchmark"
    )
    benchmark_sub = benchmark.add_subparsers(dest="benchmark", required=True)
    benchmark_sub.add_parser(
        "poiseuille", help="run the source04 exact parallel-plate benchmark"
    )

    phase3 = sub.add_parser("phase3", help="run locked Phase-3 pre-pair operations")
    phase3_sub = phase3.add_subparsers(dest="phase3_operation", required=True)
    phase3_sub.add_parser(
        "audit",
        help="run the complete 504-state numerical uncertainty audit",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "info":
        print(f"SHIFT {PROJECT_VERSION}")
        print("Shear-equivalent Hemodynamic Invariance Falsification Test")
        print(CONTRACT_SUMMARY)
        return 0

    if args.command == "doctor":
        print(f"SHIFT {PROJECT_VERSION}")
        print(f"Python: {platform.python_version()}")
        print(f"Implementation: {platform.python_implementation()}")
        print(f"Platform: {platform.platform()}")
        return 0

    if args.command == "benchmark" and args.benchmark == "poiseuille":
        from .benchmarks import representative_poiseuille_pair

        result = representative_poiseuille_pair()
        print(json.dumps(result.as_dict(), indent=2, sort_keys=True))
        return 0

    if args.command == "phase3" and args.phase3_operation == "audit":
        from .phase3_audit import run_phase3_uncertainty_audit

        result = run_phase3_uncertainty_audit()
        print(json.dumps(result.as_dict(), indent=2, sort_keys=True))
        return 0 if result.gate_passed else 2

    raise AssertionError("unreachable")


if __name__ == "__main__":
    sys.exit(main())
