#!/usr/bin/env python3
"""Regression for the q64 same-side whole-link insertion."""

from __future__ import annotations

from itertools import combinations
from math import comb, sqrt
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent
    / "open_problem_forr4_passive_floor_consolidation_round_3"
    / "searches"
)
sys.path.insert(0, str(ROOT / "searches"))
sys.path.insert(0, str(ROUND3_SEARCHES))

from opposite_endpoint_orbit_scan import endpoint_moment  # noqa: E402
from q64_paper_target_gate import (  # noqa: E402
    RESERVE_TARGET,
    THRESHOLD,
    optimize,
)
from q64_same_side_whole_link_insertion import (  # noqa: E402
    artifact_text,
    coefficient_map,
    diagnostic,
    entry_coefficient,
    inserted_coefficients,
    remaining_entries,
    same_side_whole_links,
    subset_disjointness_factor,
    target_entries,
)


def mask_factorization_check(dimension: int, left: int, right: int) -> None:
    """Check the inclusion direct-sum identity and optimized norm product."""

    universe = tuple(range(dimension))
    rows = tuple(combinations(universe, left))
    columns = tuple(combinations(universe, right))
    features = tuple(
        subset
        for level in range(left + 1)
        for subset in combinations(universe, level)
    )
    feature_index = {feature: index for index, feature in enumerate(features)}
    scales = tuple(
        sqrt(comb(right, level) / comb(left, level))
        for level in range(left + 1)
    )
    row_norm = 0.0
    column_norm = 0.0
    for row_set in rows:
        row = np.zeros(len(features))
        for level, scale in enumerate(scales):
            for subset in combinations(row_set, level):
                row[feature_index[subset]] = sqrt(scale)
        row_norm = max(row_norm, float(np.linalg.norm(row)))
        for column_set in columns:
            column = np.zeros(len(features))
            for level, scale in enumerate(scales):
                for subset in combinations(column_set, level):
                    column[feature_index[subset]] = (
                        (-1) ** level / sqrt(scale)
                    )
            expected = float(set(row_set).isdisjoint(column_set))
            if not np.isclose(np.dot(row, column), expected, atol=2e-14):
                raise AssertionError(
                    ("subset disjointness identity", left, right)
                )
            column_norm = max(column_norm, float(np.linalg.norm(column)))
    factor = subset_disjointness_factor(left, right)
    if not np.isclose(row_norm * column_norm, factor, atol=2e-13):
        raise AssertionError(("subset disjointness factor", left, right))


def endpoint_magnitude_check(order: int) -> int:
    """Check the whole 1--3 link magnitude used by the scalar reduction."""

    dimension = order * order
    checked = 0
    for support in combinations(range(dimension), 3):
        for singleton in range(0, dimension, max(1, dimension // 7)):
            moment = endpoint_moment(support, singleton, order, 3, False)
            if abs(moment) > 1 / order + 2e-14:
                raise AssertionError(
                    ("whole singleton-cubic link", order, support, moment)
                )
            checked += 1
        if checked >= 10_000:
            break
    return checked


def main() -> None:
    result = diagnostic()
    entries = target_entries()
    remaining = remaining_entries()
    if len(entries) != 96 or len(remaining) != 80:
        raise AssertionError(("same-side whole-link partition", len(entries)))
    if set(entries).intersection(remaining):
        raise AssertionError("same-side whole-link overlap")
    for entry in entries:
        links = same_side_whole_links(entry)
        if not links or not set(links).issubset({(1, 1), (1, 3), (3, 1)}):
            raise AssertionError(("unsupported same-side link", entry, links))
        if not np.isclose(coefficient_map()[entry], entry_coefficient(entry)):
            raise AssertionError(("entry coefficient", entry))
    if any(same_side_whole_links(entry) for entry in remaining):
        raise AssertionError("remaining entry has a same-side whole link")

    for degree, depth in ((5, 1), (5, 2), (7, 1), (7, 2), (7, 3)):
        mask_factorization_check(degree, depth, degree - depth)
    endpoint_checks = sum(endpoint_magnitude_check(order) for order in (4, 8))

    discrete = (
        result.closed_entries,
        result.closed_orbits,
        result.singleton_singleton_entries,
        result.singleton_cubic_entries,
        result.degree_five_extreme_entries,
        result.degree_five_balanced_entries,
        result.degree_seven_extreme_entries,
        result.degree_seven_two_five_entries,
        result.degree_seven_three_four_entries,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
    )
    if discrete != (96, 24, 28, 68, 8, 28, 24, 24, 12, 712, 808, 80):
        raise AssertionError(("same-side whole-link discrete", discrete))
    observed = (
        result.degree_five_extreme_mask_factor,
        result.degree_five_balanced_mask_factor,
        result.degree_seven_extreme_mask_factor,
        result.degree_seven_two_five_mask_factor,
        result.degree_seven_three_four_mask_factor,
        result.link_factor,
        result.minimum_coefficient,
        result.maximum_coefficient,
        result.previous_routing.total,
        result.inserted_routing.total,
        result.inserted_routing.margin_to_one_third,
        result.routing_improvement,
        result.reserve_after_declared_allowance,
        result.adaptive_multiplier_cap_retaining_allowance,
        result.remaining_class_reserve_gate,
    )
    expected = (
        3.0,
        1 + sqrt(6) + sqrt(3),
        1 + sqrt(6),
        1 + sqrt(10) + sqrt(10),
        1 + sqrt(12) + sqrt(18) + 2,
        0.015625,
        0.046875,
        0.167292848473,
        0.3284774211661729,
        0.32336258287130804,
        0.00997075046202528,
        0.005114838294864865,
        0.00897075046202528,
        (THRESHOLD - RESERVE_TARGET) / 0.32336258287130804,
        0.190775718804,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("same-side whole-link numeric", observed))
    trial = inserted_coefficients()
    for entry in remaining:
        trial[entry] = result.remaining_class_reserve_gate
    if not np.isclose(
        optimize(mapped_coefficients=trial).total,
        THRESHOLD - RESERVE_TARGET,
        rtol=0,
        atol=4e-10,
    ):
        raise AssertionError("remaining reserve gate identity")

    committed = (
        ROOT / "artifacts" / "q64_same_side_whole_link_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale same-side whole-link artifact")
    print(
        "q64 same-side whole-link insertion passed: "
        f"entries={result.closed_entries},"
        f"links={result.singleton_singleton_entries}/"
        f"{result.singleton_cubic_entries},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.inserted_routing.total:.12g},"
        f"margin={result.inserted_routing.margin_to_one_third:.12g},"
        f"remaining_open={result.remaining_open_entries},"
        f"remaining_gate={result.remaining_class_reserve_gate:.12g},"
        f"endpoint_checks={endpoint_checks},"
        "status=local_96_entry_theorem_cumulative_values_withdrawn"
    )


if __name__ == "__main__":
    main()
