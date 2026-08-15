#!/usr/bin/env python3
"""Regression for the corrected general-theory exponent ledger."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from general_theory_exponent_audit import (
    LEVELS,
    accepted_transcript_floor,
    audit,
    corrected_uniform_floor,
    one_eighth_floor,
    second_suppression_floor,
)


def main() -> None:
    witnesses, ledgers = audit()
    if tuple(witness.vertices for witness in witnesses) != LEVELS:
        raise AssertionError("witness level inventory changed")
    if any(
        witness.projective_sigma != 1 or witness.assigned_sigma != 1
        for witness in witnesses
    ):
        raise AssertionError("current interface saturator disappeared")

    expected_extras = tuple(
        Fraction(0) if level <= 8 else Fraction(level - 8, 16)
        for level in LEVELS
    )
    actual_extras = tuple(
        ledger.n_one_sixteenth_extra for ledger in ledgers
    )
    if actual_extras != expected_extras:
        raise AssertionError(("one-sixteenth level ledger", actual_extras))

    if accepted_transcript_floor() != Fraction(1, 24):
        raise AssertionError("accepted floor")
    if corrected_uniform_floor(Fraction(1, 8), 2) != Fraction(1, 16):
        raise AssertionError("stale duplicate-charge rung")
    if corrected_uniform_floor(Fraction(1, 8), 1) != Fraction(1, 8):
        raise AssertionError("corrected uniform rung")
    if second_suppression_floor() != Fraction(1, 16):
        raise AssertionError("second-suppression floor")
    if one_eighth_floor() != Fraction(1, 8):
        raise AssertionError("one-eighth floor")

    print(
        "general-theory exponent audit passed: "
        f"levels={LEVELS[0]}-{LEVELS[-1]},"
        f"accepted={accepted_transcript_floor()},"
        "stale_uniform=1/16,corrected_uniform=1/8,"
        f"second_suppression={second_suppression_floor()},"
        "interface_sigma=1"
    )


if __name__ == "__main__":
    main()
