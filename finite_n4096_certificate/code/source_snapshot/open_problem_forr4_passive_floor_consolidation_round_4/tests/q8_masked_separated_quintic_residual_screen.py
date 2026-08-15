#!/usr/bin/env python3
"""Regression for the targeted q=8 separated-quintic residual screen."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from masked_translation_subspace_screen import artifact_text  # noqa: E402
from q8_masked_separated_quintic_residual_screen import (  # noqa: E402
    residual_screen,
)


EXPECTED = {
    (1, 1, 0, 4): 6.7250649e-05,
    (2, 0, 1, 3): 7.53097278296e-05,
    (2, 1, 0, 3): 6.0607366e-05,
}


def main() -> None:
    result = residual_screen()
    if (
        result.order,
        result.templates_available,
        result.templates_screened,
        result.subspace_dimension,
        result.trials_per_template,
        result.templates_above_one,
    ) != (8, 97, 3, 6, 10, 0):
        raise AssertionError("q8 separated-quintic residual inventory")
    rows = {row.split: row for row in result.rows}
    if set(rows) != set(EXPECTED):
        raise AssertionError(("unexpected residual splits", set(rows)))
    for split, expected in EXPECTED.items():
        row = rows[split]
        if row.profile != (5, 1, 1, 5):
            raise AssertionError(("unexpected residual profile", row))
        if not np.isclose(
            row.optimized_nuclear_lower,
            expected,
            rtol=0,
            atol=6e-13,
        ):
            raise AssertionError((split, row.optimized_nuclear_lower))
        if abs(row.tangent_upper_on_orbit - row.optimized_nuclear_lower) > 2e-11:
            raise AssertionError(("translation-orbit optimization gap", split))
    committed = (
        ROOT
        / "artifacts"
        / "q8_masked_separated_quintic_residual_screen.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q8 separated-quintic residual artifact")
    print(
        "q8 separated-quintic residual regression passed: "
        f"templates={result.templates_screened},"
        f"maximum_lower={result.maximum_optimized_lower:.12g},"
        f"above_one={result.templates_above_one},"
        "status=no_counterexample_diagnostic_only"
    )


if __name__ == "__main__":
    main()
