#!/usr/bin/env python3
"""Regression for selected mixed translation-orbit shape screens."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from masked_translation_mixture_screen import (  # noqa: E402
    artifact_text,
    diagnostic,
)


INDICES = (7, 18, 20, 23, 36, 42, 46, 49)


def main() -> None:
    q4 = diagnostic(
        order=4,
        orbit_shapes=3,
        group_dimension=5,
        selected_indices=INDICES,
    )
    q8 = diagnostic(
        order=8,
        orbit_shapes=3,
        group_dimension=5,
        selected_indices=INDICES,
    )
    if (q4.templates_screened, q8.templates_screened) != (8, 8):
        raise AssertionError("mixed translation selected inventories")
    if q4.templates_above_one or q8.templates_above_one:
        raise AssertionError("unexpected mixed translation counterexample")
    expected = (0.019404557550692126, 0.0002790178544551007)
    observed = (q4.maximum_optimized_lower, q8.maximum_optimized_lower)
    if not np.allclose(observed, expected, rtol=0, atol=3e-14):
        raise AssertionError(("mixed translation maxima", observed))
    if q4.maximum_shape_tangent_upper - q4.maximum_optimized_lower > 4e-5:
        raise AssertionError("q4 mixed-shape optimization gap")
    if q8.maximum_shape_tangent_upper - q8.maximum_optimized_lower > 4e-12:
        raise AssertionError("q8 mixed-shape optimization gap")
    for order, result in ((4, q4), (8, q8)):
        committed = (
            ROOT
            / "artifacts"
            / f"q{order}_masked_translation_mixture_screen.json"
        ).read_text(encoding="utf-8")
        if committed != artifact_text(result):
            raise AssertionError(("stale mixed translation artifact", order))
    print(
        "masked translation-mixture regressions passed: "
        f"q4_max={q4.maximum_optimized_lower:.12g},"
        f"q4_gap={q4.maximum_shape_tangent_upper - q4.maximum_optimized_lower:.12g},"
        f"q8_max={q8.maximum_optimized_lower:.12g},"
        f"q8_gap={q8.maximum_shape_tangent_upper - q8.maximum_optimized_lower:.12g},"
        "status=shape_mixing_visible_no_counterexample"
    )


if __name__ == "__main__":
    main()
