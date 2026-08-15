#!/usr/bin/env python3
"""Regression for the all-template q=4 and selected q=8 masked screens."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from masked_translation_subspace_screen import (  # noqa: E402
    artifact_text,
    diagnostic,
)


Q8_INDICES = (7, 18, 20, 23, 36, 42, 46, 49)


def main() -> None:
    q4 = diagnostic(order=4, subspace_dimension=6, trials=3)
    if (
        q4.templates_available,
        q4.templates_screened,
        q4.templates_above_one,
    ) != (97, 97, 0):
        raise AssertionError("q4 masked translation inventory")
    if not np.isclose(
        q4.maximum_optimized_lower,
        0.02209708691207961,
        rtol=0,
        atol=3e-14,
    ):
        raise AssertionError(("q4 maximum masked lower", q4.maximum_optimized_lower))
    if max(
        abs(row.tangent_upper_on_orbit - row.optimized_nuclear_lower)
        for row in q4.rows
    ) > 2e-10:
        raise AssertionError("q4 orbit laws are not at the concavity tangent")
    q4_riskiest = max(q4.rows, key=lambda row: row.optimized_nuclear_lower)
    if (q4_riskiest.profile, q4_riskiest.split) != (
        (1, 1, 9, 1),
        (0, 1, 5, 0),
    ):
        raise AssertionError(("q4 riskiest template", q4_riskiest))

    q8 = diagnostic(
        order=8,
        subspace_dimension=6,
        trials=3,
        selected_indices=Q8_INDICES,
    )
    if (q8.templates_screened, q8.templates_above_one) != (8, 0):
        raise AssertionError("q8 selected masked translation inventory")
    if not np.isclose(
        q8.maximum_optimized_lower,
        6.01796939448103e-05,
        rtol=0,
        atol=3e-16,
    ):
        raise AssertionError(("q8 maximum masked lower", q8.maximum_optimized_lower))
    if max(
        abs(row.tangent_upper_on_orbit - row.optimized_nuclear_lower)
        for row in q8.rows
    ) > 2e-11:
        raise AssertionError("q8 orbit laws are not at the concavity tangent")

    committed_q4 = (
        ROOT / "artifacts" / "q4_masked_translation_subspace_screen.json"
    ).read_text(encoding="utf-8")
    if committed_q4 != artifact_text(q4):
        raise AssertionError("stale q4 masked translation artifact")
    committed_q8 = (
        ROOT / "artifacts" / "q8_masked_translation_subspace_screen.json"
    ).read_text(encoding="utf-8")
    if committed_q8 != artifact_text(q8):
        raise AssertionError("stale q8 masked translation artifact")
    print(
        "masked translation-subspace regressions passed: "
        f"q4_templates={q4.templates_screened},"
        f"q4_max={q4.maximum_optimized_lower:.12g},"
        f"q8_templates={q8.templates_screened},"
        f"q8_max={q8.maximum_optimized_lower:.12g},"
        "status=no_counterexample_translation_twirling_mechanism_identified"
    )


if __name__ == "__main__":
    main()
