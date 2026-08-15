#!/usr/bin/env python3
"""Exact dose-six constants for unordered joint occurrence packing."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb, sqrt
import random


SEED = 2026071417
BLOCKS = 4
DOSE = 6
TWO_COPY_DOSE = 2 * DOSE
BETA = Fraction(5, 6)


def profiles() -> list[tuple[int, ...]]:
    return [
        profile
        for profile in product(range(1, TWO_COPY_DOSE + 1, 2), repeat=BLOCKS)
        if sum(profile) <= TWO_COPY_DOSE
    ]


def capacity_states() -> list[tuple[int, ...]]:
    return [
        state
        for state in product(range(TWO_COPY_DOSE + 1), repeat=BLOCKS)
        if sum(state) <= TWO_COPY_DOSE
    ]


def packing_constant(
    profile: tuple[int, ...], states: list[tuple[int, ...]]
) -> tuple[int, tuple[int, ...]]:
    return max(
        (
            comb(state[0], profile[0])
            * comb(state[1], profile[1])
            * comb(state[2], profile[2])
            * comb(state[3], profile[3]),
            state,
        )
        for state in states
    )


def elementary_coefficient(multiplicities: list[int], degree: int) -> int:
    coefficients = [0 for _ in range(degree + 1)]
    coefficients[0] = 1
    for multiplicity in multiplicities:
        for current in range(degree, 0, -1):
            coefficients[current] += multiplicity * coefficients[current - 1]
    return coefficients[degree]


def random_incidence_checks(rng: random.Random) -> int:
    checked = 0
    for universe_size in range(1, 10):
        for _ in range(80):
            container_count = rng.randrange(1, 9)
            multiplicities = [0 for _ in range(universe_size)]
            for _ in range(container_count):
                support = [
                    coordinate
                    for coordinate in range(universe_size)
                    if rng.randrange(2)
                ]
                for coordinate in support:
                    multiplicities[coordinate] += 1
            total_slots = sum(multiplicities)
            for degree in range(universe_size + 1):
                direct = elementary_coefficient(multiplicities, degree)
                bound = comb(total_slots, degree)
                if direct > bound:
                    raise AssertionError(
                        ("unordered slot packing", multiplicities, degree)
                    )
                checked += 1
    return checked


def main() -> None:
    rng = random.Random(SEED)
    incidence_checks = random_incidence_checks(rng)
    values = profiles()
    states = capacity_states()
    if len(values) != 70:
        raise AssertionError(("profile count", len(values)))

    rows_by_degree: dict[int, list[tuple[int, tuple[int, ...], tuple[int, ...]]]] = {}
    attenuated_sums: dict[int, float] = {}
    for profile in values:
        value, state = packing_constant(profile, states)
        degree = sum(profile)
        rows_by_degree.setdefault(degree, []).append((value, profile, state))
        attenuated_sums[degree] = attenuated_sums.get(degree, 0.0) + (
            sqrt(value) * float(BETA**degree)
        )

    expected_profile_counts = {4: 1, 6: 4, 8: 10, 10: 20, 12: 35}
    expected_maximum_mass = {4: 81, 6: 160, 8: 126, 10: 36, 12: 1}
    if {
        degree: len(rows) for degree, rows in rows_by_degree.items()
    } != expected_profile_counts:
        raise AssertionError(("degree profile counts", rows_by_degree))
    for degree, rows in rows_by_degree.items():
        maximum = max(row[0] for row in rows)
        if maximum != expected_maximum_mass[degree]:
            raise AssertionError(("packing maximum", degree, maximum))

    print(
        "joint occurrence profile packing passed: "
        f"incidence_checks={incidence_checks},profiles={len(values)},"
        f"capacity_states={len(states)}"
    )
    for degree in sorted(rows_by_degree):
        rows = rows_by_degree[degree]
        maximum = max(rows)
        print(
            f"degree={degree},profiles={len(rows)},"
            f"max_square_mass={maximum[0]},"
            f"max_sqrt_mass={sqrt(maximum[0]):.12g},"
            f"maximizing_profile={maximum[1]},"
            f"capacity={maximum[2]},"
            f"beta_weighted_sqrt_sum={attenuated_sums[degree]:.12g}"
        )


if __name__ == "__main__":
    main()
