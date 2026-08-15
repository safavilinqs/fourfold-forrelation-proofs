#!/usr/bin/env python3
"""Regression for exact q=4/q=8 masked translation covariance."""

from __future__ import annotations

from pathlib import Path
import random
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from masked_translation_covariance import (  # noqa: E402
    chain_translation_sign,
    disjoint_union,
    link_translation_sign,
    separated_occurrence_signs,
    translation_cocycle,
    translate_configuration,
    translate_support,
    xor_shifts,
)
from signed_permutation_link_moment import (  # noqa: E402
    chain_moment,
    link_moment,
)


PROFILES = (
    (3, 3, 3, 3),
    (1, 1, 3, 7),
    (1, 1, 5, 5),
    (1, 1, 9, 1),
    (1, 3, 3, 5),
)


def random_support(
    order: int,
    degree: int,
    rng: random.Random,
) -> tuple[int, ...]:
    return tuple(sorted(rng.sample(range(order * order), degree)))


def check_link(order: int, rng: random.Random, trials: int) -> int:
    comparisons = 0
    for _ in range(trials):
        left = random_support(order, rng.choice((1, 3, 5, 7)), rng)
        right = random_support(order, rng.choice((1, 3, 5, 7)), rng)
        left_shift = (rng.randrange(order), rng.randrange(order))
        right_shift = (rng.randrange(order), rng.randrange(order))
        original = link_moment(order, left, right)
        translated = link_moment(
            order,
            translate_support(order, left, left_shift),
            translate_support(order, right, right_shift),
        )
        expected = original * link_translation_sign(
            order, left, right, left_shift, right_shift
        )
        if translated != expected:
            raise AssertionError(
                ("link translation covariance", order, translated, expected)
            )
        comparisons += 1
    return comparisons


def check_chain_and_separation(
    order: int,
    rng: random.Random,
    trials: int,
) -> int:
    comparisons = 0
    for _ in range(trials):
        profile = rng.choice(PROFILES)
        supports = tuple(
            random_support(order, degree, rng) for degree in profile
        )
        shifts = tuple(
            (rng.randrange(order), rng.randrange(order)) for _ in range(4)
        )
        original = chain_moment(order, supports)
        translated = chain_moment(
            order, translate_configuration(order, supports, shifts)
        )
        expected = original * chain_translation_sign(
            order, supports, shifts
        )
        if translated != expected:
            raise AssertionError(
                ("chain translation covariance", order, translated, expected)
            )

        row = []
        column = []
        for support in supports:
            selected = rng.randrange(len(support) + 1)
            chosen = set(rng.sample(support, selected))
            row.append(tuple(value for value in support if value in chosen))
            column.append(tuple(value for value in support if value not in chosen))
        row_configuration = tuple(row)
        column_configuration = tuple(column)
        full = disjoint_union(row_configuration, column_configuration)
        row_sign, column_sign, full_sign = separated_occurrence_signs(
            order,
            row_configuration,
            column_configuration,
            shifts,
        )
        if full != supports or row_sign * column_sign != full_sign:
            raise AssertionError("occurrence translation sign separation")
        comparisons += 1
    return comparisons


def check_projective_cocycle(
    order: int,
    rng: random.Random,
    trials: int,
) -> int:
    comparisons = 0
    for _ in range(trials):
        profile = rng.choice(PROFILES)
        sizes = tuple(rng.randrange(degree + 1) for degree in profile)
        support = tuple(
            random_support(order, size, rng) for size in sizes
        )
        left = tuple(
            (rng.randrange(order), rng.randrange(order)) for _ in range(4)
        )
        right = tuple(
            (rng.randrange(order), rng.randrange(order)) for _ in range(4)
        )
        observed = (
            chain_translation_sign(order, support, right)
            * chain_translation_sign(
                order,
                translate_configuration(order, support, right),
                left,
            )
            * chain_translation_sign(
                order, support, xor_shifts(left, right)
            )
        )
        expected = translation_cocycle(
            order, sizes, left, right
        )
        if observed != expected:
            raise AssertionError(
                ("support-dependent translation cocycle", order, sizes)
            )
        comparisons += 1
    return comparisons


def main() -> None:
    rng = random.Random(2026071604)
    counts = []
    for order in (4, 8):
        counts.append(
            (
                order,
                check_link(order, rng, 400),
                check_chain_and_separation(order, rng, 400),
                check_projective_cocycle(order, rng, 400),
            )
        )
    print(
        "masked translation covariance passed: "
        + ",".join(
            f"q{order}_link={link},q{order}_chain={chain},"
            f"q{order}_cocycle={cocycle}"
            for order, link, chain, cocycle in counts
        )
        + ",status=exact_character_covariance_and_occurrence_separation"
    )


if __name__ == "__main__":
    main()
