#!/usr/bin/env python3
"""Regression for the non-invariant opposite-endpoint two-axis screen."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ROUND3 = ROOT.parent / "open_problem_forr4_passive_floor_consolidation_round_3"
sys.path.insert(0, str(ROOT / "searches"))
sys.path.insert(0, str(ROUND3 / "searches"))

from opposite_endpoint_mixed_orbit_q4 import (  # noqa: E402
    build_data as build_q4_data,
    evaluate as evaluate_q4,
)
from opposite_endpoint_orbit_scan import (  # noqa: E402
    triple_orbit_representatives,
)
from opposite_endpoint_two_axis_mixture import (  # noqa: E402
    artifact_text,
    coefficient,
    search,
)


def unreduced_q4(equal_mass: float, horizontal_mass: float) -> float:
    order = 4
    data = build_q4_data(order)
    vertical = tuple(row * order for row in range(1, order))
    horizontal = tuple(range(1, order))
    row = np.zeros((len(data.differences), len(data.differences)))
    for left in vertical:
        left_index = int(np.where(data.differences == left)[0][0])
        for right in vertical:
            right_index = int(np.where(data.differences == right)[0][0])
            probability = (
                equal_mass / (order - 1)
                if left == right
                else (1 - equal_mass) / ((order - 1) * (order - 2))
            )
            row[left_index, right_index] = probability * (1 - horizontal_mass)
        for right in horizontal:
            right_index = int(np.where(data.differences == right)[0][0])
            row[left_index, right_index] = horizontal_mass / (order - 1) ** 2

    vertical_triples = tuple(
        tuple(row_label * order for row_label in triple)
        for triple in triple_orbit_representatives(order)
    )
    column = np.zeros(len(data.triples))
    for triple in vertical_triples:
        column[data.triples.index(triple)] = 1 / len(vertical_triples)
    return evaluate_q4(data, row, column).objective


def main() -> None:
    equal_mass = 2 / 5
    horizontal_mass = 1 / 20
    reduced = coefficient(4, equal_mass, horizontal_mass)
    direct = unreduced_q4(equal_mass, horizontal_mass)
    if not np.isclose(reduced, direct, atol=4e-13):
        raise AssertionError(("q4 two-axis reduction", reduced, direct))

    result = search()
    observed = (
        result.invariant_equal_mass,
        result.invariant_coefficient,
        result.best_horizontal_mass,
        result.best_coefficient,
        result.gate_headroom,
    )
    expected = (
        0.0192515087131167,
        0.0395996495753725,
        0.000647887695545267,
        0.0396118487000776,
        0.00185046959643697,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=4e-12):
        raise AssertionError(("q32 two-axis search", observed))
    if result.best_coefficient >= result.physical_gate:
        raise AssertionError(("two-axis family crossed gate", result))
    if result.decision != "does_not_cross_leading_orbit_kill_gate":
        raise AssertionError(("two-axis decision", result.decision))

    committed = (
        ROOT / "artifacts" / "opposite_endpoint_two_axis_screen.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale opposite-endpoint two-axis artifact")

    print(
        "opposite-endpoint two-axis screen passed: "
        f"q4_direct={direct:.12g},"
        f"invariant={result.invariant_coefficient:.12g},"
        f"horizontal_mass={result.best_horizontal_mass:.12g},"
        f"best={result.best_coefficient:.12g},"
        f"gate={result.physical_gate:.12g},"
        f"headroom={result.gate_headroom:.12g},"
        f"decision={result.decision}"
    )


if __name__ == "__main__":
    main()
