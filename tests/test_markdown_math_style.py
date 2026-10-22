"""Regression tests for the public Markdown math convention."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_MARKDOWN_MATH_DELIMITERS = (r"\(", r"\)", r"\[", r"\]")


def test_markdown_math_uses_dollar_delimiters_only() -> None:
    offenders: list[str] = []

    for path in sorted(ROOT.rglob("*.md")):
        text = path.read_text(encoding="utf-8")
        relative = path.relative_to(ROOT)
        for delimiter in FORBIDDEN_MARKDOWN_MATH_DELIMITERS:
            if delimiter in text:
                offenders.append(f"{relative}: contains forbidden {delimiter!r}")

    assert not offenders, (
        "SHIFT Markdown math must use $...$ inline and $$...$$ for display math.\n"
        + "\n".join(offenders)
    )
