#!/usr/bin/env python3
"""Regression for the explicit N=4096 active six-dose resource row."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from active_six_resource_row import (  # noqa: E402
    artifact_text,
    resource_row,
)


def main() -> None:
    result = resource_row()
    observed = (
        result.dimension,
        result.sign_blocks,
        result.sign_modes,
        result.flags,
        result.single_photons,
        result.logical_path_dimension_per_flag,
        result.logical_mode_dimension_per_flag,
        result.hard_dose_per_flag,
        result.total_hard_dose,
        result.left_branch.charged_sign_blocks,
        result.left_branch.public_hadamards,
        result.right_branch.charged_sign_blocks,
        result.right_branch.public_hadamards,
        result.postselection,
    )
    expected = (
        4096,
        4,
        16_384,
        3,
        3,
        2,
        4096,
        2,
        6,
        (1, 2),
        1,
        (4, 3),
        2,
        False,
    )
    if observed != expected:
        raise AssertionError(("active resource ledger", observed))
    if result.left_branch.hard_dose != 2 or result.right_branch.hard_dose != 2:
        raise AssertionError("active branch hard-dose mismatch")
    if Fraction(result.flag_correct_probability_at_promise_boundary) != Fraction(5, 8):
        raise AssertionError("active flag endpoint probability")
    if Fraction(result.majority_error_exact) != Fraction(81, 256):
        raise AssertionError("active majority error")
    if Fraction(result.margin_below_one_third_exact) != Fraction(13, 768):
        raise AssertionError("active error margin")
    if not result.majority_error < 1 / 3:
        raise AssertionError("active protocol misses the decision threshold")

    committed = (ROOT / "artifacts" / "active_six_resource_row.json").read_text(
        encoding="utf-8"
    )
    if committed != artifact_text(result):
        raise AssertionError("stale active six-dose resource artifact")

    print(
        "active six-dose resource row passed: "
        f"N={result.dimension},"
        f"M={result.sign_modes},"
        f"flags={result.flags},"
        f"dose={result.total_hard_dose},"
        f"error={result.majority_error_exact},"
        f"margin={result.margin_below_one_third_exact}"
    )


if __name__ == "__main__":
    main()
