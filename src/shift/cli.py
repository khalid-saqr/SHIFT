"""Public command-line entry point for SHIFT."""

import argparse
import platform
import sys
from importlib.metadata import version

PROJECT_VERSION = version("shift-hemodynamics")
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

    raise AssertionError("unreachable")


if __name__ == "__main__":
    sys.exit(main())
