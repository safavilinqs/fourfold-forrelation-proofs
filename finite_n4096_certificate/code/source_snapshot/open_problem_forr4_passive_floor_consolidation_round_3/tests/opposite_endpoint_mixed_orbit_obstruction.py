#!/usr/bin/env python3
"""Regression for the opposite-endpoint mixed-orbit physical witness."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from opposite_endpoint_mixed_orbit_q4 import build_data, evaluate
from opposite_endpoint_vertical_mixture_witness import (
    forced_ledger_euclidean_repair,
    forced_ledger_obstruction,
    forced_split_coefficients,
    vertical_mixture_witness,
)


def reduced_law_check(order: int, expected: float) -> None:
    data = build_data(order)
    vertical = np.flatnonzero(data.differences % order == 0)
    row = np.zeros((len(data.differences), len(data.differences)))
    row[np.ix_(vertical, vertical)] = 1 / len(vertical) ** 2
    same_column = tuple(
        index
        for index, triple in enumerate(data.triples)
        if len({coordinate % order for coordinate in triple}) == 1
    )
    expected_orbits = (order - 1) * (order - 2) // 6
    if len(same_column) != expected_orbits:
        raise AssertionError(
            ("same-column triple orbits", order, len(same_column))
        )
    column = np.zeros(len(data.triples))
    column[list(same_column)] = 1 / len(same_column)
    observed = evaluate(data, row, column).objective
    if not np.isclose(observed, expected, atol=2e-13):
        raise AssertionError(
            ("mixed-orbit reduction", order, observed, expected)
        )


def main() -> None:
    expected = {
        4: 0.0676582346706592,
        8: 0.127869435555387,
        16: 0.0853097508344982,
        32: 0.039593955294628,
    }
    witnesses = {
        order: vertical_mixture_witness(order) for order in expected
    }
    for order, target in expected.items():
        observed = witnesses[order].coefficient
        if not np.isclose(observed, target, atol=2e-13):
            raise AssertionError(
                ("vertical mixture coefficient", order, observed, target)
            )
    reduced_law_check(4, expected[4])
    reduced_law_check(8, expected[8])

    coefficients = forced_split_coefficients(expected[32])
    if len(coefficients) != 4 or any(
        not np.isclose(value, expected[32])
        for value in coefficients.values()
    ):
        raise AssertionError(("forced critical cuts", coefficients))
    obstruction = forced_ledger_obstruction(expected[32])
    if not 0.78246 < obstruction.optimal_beta < 0.78248:
        raise AssertionError(("post-repair forced-cut beta", obstruction))
    if not 0.28378 < obstruction.optimized_total < 0.28379:
        raise AssertionError(("post-repair forced-cut total", obstruction))
    if not obstruction.threshold_overshoot < -0.0495:
        raise AssertionError(("post-repair forced-cut slack", obstruction))
    repair = forced_ledger_euclidean_repair(expected[32])
    if not 0.78114 < repair.optimal_beta < 0.78116:
        raise AssertionError(("repaired forced-cut beta", repair))
    if not 0.28262 < repair.optimized_total < 0.28264:
        raise AssertionError(("repaired forced-cut total", repair))
    if not repair.threshold_slack > 0.0507:
        raise AssertionError(("repaired forced-cut slack", repair))
    print(
        "opposite endpoint mixed-orbit witness passed: "
        f"q32_coefficient={expected[32]:.12g},"
        f"optimal_beta={obstruction.optimal_beta:.12g},"
        f"optimized_total={obstruction.optimized_total:.12g},"
        f"post_repair_slack={-obstruction.threshold_overshoot:.12g},"
        f"repaired_beta={repair.optimal_beta:.12g},"
        f"repaired_total={repair.optimized_total:.12g},"
        f"repaired_slack={repair.threshold_slack:.12g}"
    )


if __name__ == "__main__":
    main()
