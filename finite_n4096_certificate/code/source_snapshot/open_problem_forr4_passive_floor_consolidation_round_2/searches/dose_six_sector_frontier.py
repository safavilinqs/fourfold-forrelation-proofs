#!/usr/bin/env python3
"""Enumerate the exact odd Fourier-profile frontier at hard dose six."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb


DOSE = 6
BETA = Fraction(5, 6)


def profiles() -> list[tuple[int, int, int, int]]:
    result = []
    for excess in product(range(DOSE), repeat=4):
        if sum(excess) <= DOSE - 2:
            result.append(tuple(1 + 2 * value for value in excess))
    return result


def record_sector_count(profile: tuple[int, int, int, int]) -> int:
    result = 1
    for left, right in zip(profile, profile[1:]):
        maximum = min(left, right)
        result *= (maximum + 1) // 2
    return result


def disposition(profile: tuple[int, int, int, int]) -> str:
    decorated = [index for index, degree in enumerate(profile) if degree > 1]
    if not decorated:
        return "minimal_closed"
    if len(decorated) == 1:
        block = decorated[0]
        if block in (0, 3):
            return "endpoint_closed"
        if profile[block] == 3:
            return "middle_cubic_local_closed"
        if profile[block] == 5:
            return "middle_quintic_local_closed"
        return "middle_higher_open"
    if sum(profile) == 8 and len(decorated) == 2:
        pair = tuple(decorated)
        if pair in ((0, 2), (1, 2), (1, 3)):
            return "double_cubic_entry_bounded"
        if pair in ((0, 1), (0, 3), (2, 3)):
            return "double_cubic_occurrence_bounded"
    return "multiple_decorations_open"


def main() -> None:
    values = profiles()
    if len(values) != comb(8, 4):
        raise AssertionError(("odd profile count", len(values)))

    by_degree: dict[int, dict[str, int]] = {}
    weighted: dict[str, Fraction] = {}
    total_record_sectors = 0
    for profile in values:
        degree = sum(profile)
        state = disposition(profile)
        by_degree.setdefault(degree, {})[state] = (
            by_degree.setdefault(degree, {}).get(state, 0) + 1
        )
        weighted[state] = weighted.get(state, Fraction()) + BETA**degree
        total_record_sectors += record_sector_count(profile)

    closed = {
        "minimal_closed",
        "endpoint_closed",
        "middle_cubic_local_closed",
        "middle_quintic_local_closed",
        "double_cubic_entry_bounded",
        "double_cubic_occurrence_bounded",
    }
    closed_profiles = sum(
        count
        for states in by_degree.values()
        for state, count in states.items()
        if state in closed
    )
    if closed_profiles != 19:
        raise AssertionError(("closed profile count", closed_profiles))
    if total_record_sectors != 130:
        raise AssertionError(("record triple count", total_record_sectors))
    if any(
        state not in closed
        for degree in (4, 6)
        for state in by_degree[degree]
    ):
        raise AssertionError(("degree-six local frontier", by_degree))

    rows = []
    for degree in sorted(by_degree):
        states = ",".join(
            f"{state}={count}"
            for state, count in sorted(by_degree[degree].items())
        )
        rows.append(f"degree={degree}:profiles={sum(by_degree[degree].values())},{states}")
    weighted_rows = ",".join(
        f"{state}={float(value):.12g}"
        for state, value in sorted(weighted.items())
    )
    print(
        "dose-six signed-permutation sector frontier:\n"
        + "\n".join(rows)
        + "\n"
        + f"total_profiles={len(values)},locally_closed={closed_profiles},"
        + f"open={len(values)-closed_profiles},record_triples={total_record_sectors}\n"
        + f"beta=5/6 weighted_profile_counts:{weighted_rows}"
    )


if __name__ == "__main__":
    main()
