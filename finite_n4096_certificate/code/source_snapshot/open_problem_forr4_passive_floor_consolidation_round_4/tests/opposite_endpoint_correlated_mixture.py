#!/usr/bin/env python3
"""Regression for the invariant opposite-endpoint correlation search."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ROUND3 = ROOT.parent / "open_problem_forr4_passive_floor_consolidation_round_3"
sys.path.insert(0, str(ROOT / "searches"))
sys.path.insert(0, str(ROUND3 / "searches"))

from opposite_endpoint_correlated_mixture import (  # noqa: E402
    correlated_vertical_mixture,
    search_correlated_vertical_mixture,
)
from opposite_endpoint_mixed_orbit_q4 import build_data, evaluate  # noqa: E402


def q4_reduction_check(equal_mass: float) -> None:
    """Compare one invariant law with the inherited unreduced q=4 sum."""

    order = 4
    invariant_law = np.full(
        (order - 1, order - 1),
        (1 - equal_mass) / ((order - 1) * (order - 2)),
    )
    np.fill_diagonal(invariant_law, equal_mass / (order - 1))
    reduced = correlated_vertical_mixture(order, equal_mass).coefficient

    data = build_data(order)
    vertical_indices = np.flatnonzero(data.differences % order == 0)
    full_row = np.zeros((len(data.differences), len(data.differences)))
    full_row[np.ix_(vertical_indices, vertical_indices)] = invariant_law
    same_column = tuple(
        index
        for index, triple in enumerate(data.triples)
        if len({coordinate % order for coordinate in triple}) == 1
    )
    full_column = np.zeros(len(data.triples))
    full_column[list(same_column)] = 1 / len(same_column)
    direct = evaluate(data, full_row, full_column).objective
    if not np.isclose(reduced, direct, atol=3e-13):
        raise AssertionError(("q4 invariant reduction", reduced, direct))


def main() -> None:
    q4_reduction_check(1 / 3)
    q4_reduction_check(2 / 5)
    result = search_correlated_vertical_mixture()
    expected = (
        0.039593955294628,
        0.0192515195006088,
        0.0395996495753725,
    )
    observed = (
        result.independent_coefficient,
        result.best_equal_mass,
        result.best_coefficient,
    )
    if not np.allclose(observed, expected, atol=4e-13):
        raise AssertionError(("q32 invariant correlation search", result))
    physical_gate = 0.0414623182965146
    if not result.best_coefficient < physical_gate:
        raise AssertionError(("invariant family unexpectedly obstructs gate", result))

    print(
        "opposite-endpoint invariant correlation search passed: "
        f"independent={result.independent_coefficient:.12g},"
        f"best_mass={result.best_equal_mass:.12g},"
        f"best={result.best_coefficient:.12g},"
        f"gate={physical_gate:.12g},"
        f"headroom={physical_gate-result.best_coefficient:.12g}"
    )


if __name__ == "__main__":
    main()
