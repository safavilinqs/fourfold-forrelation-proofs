#!/usr/bin/env python3
"""Regression for the open high-degree record/incidence frontier."""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from high_degree_record_incidence_frontier import (
    coarse_target_priorities,
    compatible_record_triples,
    dose_six_relevant_entries,
    frontier_summary,
    open_profiles,
    record_incidence_for_order,
    symmetry_orbits,
)


def direct_one_axis_incidence(
    order: int, degree: int, selected: int, record: int
) -> int:
    """Brute-force the small-order incidence maximum."""

    cells = tuple((row, column) for row in range(order) for column in range(order))
    incidence: dict[tuple[tuple[int, int], ...], int] = defaultdict(int)
    for support in combinations(cells, degree):
        odd_columns = sum(
            sum(cell[1] == column for cell in support) % 2
            for column in range(order)
        )
        if odd_columns != record:
            continue
        for fixed in combinations(support, selected):
            incidence[fixed] += 1
    return max(incidence.values(), default=0)


def main() -> None:
    for degree in (1, 3, 5):
        for selected in range(degree + 1):
            for record in range(1, min(4, degree) + 1, 2):
                direct = direct_one_axis_incidence(
                    4, degree, selected, record
                )
                dynamic = record_incidence_for_order(
                    4, degree, selected, record
                )
                if direct != dynamic:
                    raise AssertionError(
                        (
                            "record incidence mismatch",
                            degree,
                            selected,
                            record,
                            direct,
                            dynamic,
                        )
                    )

    profiles = open_profiles()
    if len(profiles) != 43:
        raise AssertionError(("open profile count", len(profiles)))
    if sum(len(compatible_record_triples(profile)) for profile in profiles) != 92:
        raise AssertionError("compatible record-triple inventory changed")

    summary = frontier_summary()
    expected = (
        6016,
        14624,
        1138,
        4878,
        256,
        896,
        882,
        5120,
        888,
        0,
    )
    actual = (
        summary.profile_splits,
        summary.record_sectors,
        summary.certified_profile_splits,
        summary.unresolved_profile_splits,
        summary.certified_degree_ten,
        summary.degree_ten_splits,
        summary.certified_degree_twelve,
        summary.degree_twelve_splits,
        summary.dose_six_relevant_splits,
        summary.certified_dose_six_relevant_splits,
    )
    if actual != expected:
        raise AssertionError(("frontier summary", actual))

    relevant = dose_six_relevant_entries()
    orbits = symmetry_orbits(relevant)
    if set().union(*map(set, orbits)) != set(relevant):
        raise AssertionError("symmetry orbits do not cover the relevant frontier")
    if sum(map(len, orbits)) != len(relevant):
        raise AssertionError("symmetry orbits overlap")

    priorities = coarse_target_priorities(limit=1)
    leading = priorities[0]
    if not np.isclose(
        leading.perron_contribution, 0.00394518332418, atol=3e-12
    ):
        raise AssertionError(("leading route contribution", leading))
    if len(leading.entries) != 4 or not np.isclose(
        leading.current_coefficients[0], 0.039593955294628, atol=1e-15
    ):
        raise AssertionError(("leading physical orbit", leading))

    print(
        "high-degree record/incidence frontier passed: "
        f"profiles={summary.profiles},"
        f"record_sectors={summary.record_sectors},"
        f"certified={summary.certified_profile_splits},"
        f"dose_six_relevant={summary.dose_six_relevant_splits},"
        f"dose_six_certified={summary.certified_dose_six_relevant_splits},"
        f"leading_contribution={leading.perron_contribution:.12g}"
    )


if __name__ == "__main__":
    main()
