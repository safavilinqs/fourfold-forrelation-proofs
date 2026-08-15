#!/usr/bin/env python3
"""Regression for the exact full-translation projective block screen."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from masked_translation_full_group_screen import (  # noqa: E402
    algebra_transform_check,
    artifact_text,
    cocycle_bilinearity_check,
    diagnostic,
    physical_transform_check,
)


def main() -> None:
    if algebra_transform_check() != 92:
        raise AssertionError("full-group algebra check inventory")
    if cocycle_bilinearity_check() != 4800:
        raise AssertionError("cocycle bicharacter check inventory")
    direct, irrep, characters = physical_transform_check()
    if not np.isclose(direct, 0.3535533905932742, rtol=0, atol=2e-15):
        raise AssertionError(("q2 direct masked transform", direct))
    if (irrep, characters) != (16, 1):
        raise AssertionError(("q2 projective type", irrep, characters))

    result = diagnostic()
    if (result.templates, result.coefficients_above_one) != (97, 0):
        raise AssertionError("full-group screen inventory")
    observed = (
        result.canonical_pure_maximum,
        result.focused_pure_maximum,
        result.mixed_maximum_lower,
        result.mixed_maximum_tangent_upper,
    )
    expected = (
        0.08369866694739853,
        0.17677669529663698,
        0.17677669314724448,
        0.17677669529663687,
    )
    if not np.allclose(observed, expected, rtol=0, atol=3e-14):
        raise AssertionError(("full-group maxima", observed))
    if tuple(row.normalized_rank for row in result.mixed_rows) != (0, 4, 8):
        raise AssertionError("one finite mixed simplex per projective type")
    if any(
        row.concavity_tangent_upper >= 1 for row in result.mixed_rows
    ):
        raise AssertionError("unexpected finite mixed-simplex counterexample")

    committed = (
        ROOT
        / "artifacts"
        / "q4_masked_translation_full_group_screen.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q4 full-translation artifact")
    print(
        "q4 masked full-translation regression passed: "
        f"templates={result.templates},"
        f"canonical_max={result.canonical_pure_maximum:.12g},"
        f"focused_max={result.focused_pure_maximum:.12g},"
        f"mixed_tangent={result.mixed_maximum_tangent_upper:.12g},"
        "status=exact_clifford_reduction_no_finite_counterexample"
    )


if __name__ == "__main__":
    main()
