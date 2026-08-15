#!/usr/bin/env python3
"""Regression tests for the exact sparse-support plant moment."""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations, product
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from signed_permutation_link_moment import moment  # noqa: E402


SEED = 2026071501


def sylvester(order: int) -> np.ndarray:
    result = np.asarray([[1]], dtype=np.int8)
    while len(result) < order:
        result = np.block([[result, result], [result, -result]])
    return result


def exact_pairs(order: int) -> tuple[np.ndarray, np.ndarray]:
    walsh = sylvester(order)
    left: list[np.ndarray] = []
    right: list[np.ndarray] = []
    for permutation in permutations(range(order)):
        for signs in product((-1, 1), repeat=order):
            signed = np.zeros((order, order), dtype=np.int8)
            for column, row in enumerate(permutation):
                signed[row, column] = signs[column]
            left.append((walsh @ signed).reshape(-1))
            right.append((signed @ walsh).reshape(-1))
    return np.asarray(left), np.asarray(right)


def enumerated_moment(
    left_values: np.ndarray,
    right_values: np.ndarray,
    left_support: tuple[int, ...],
    right_support: tuple[int, ...],
) -> Fraction:
    values = np.ones(len(left_values), dtype=np.int64)
    for coordinate in left_support:
        values *= left_values[:, coordinate]
    for coordinate in right_support:
        values *= right_values[:, coordinate]
    return Fraction(int(values.sum()), len(values))


def main() -> None:
    rng = np.random.default_rng(SEED)
    comparisons = 0
    for order, trials in ((2, 80), (4, 160)):
        left_values, right_values = exact_pairs(order)
        dimension = order * order
        for _ in range(trials):
            left_size = int(rng.integers(0, min(13, dimension + 1)))
            right_size = int(rng.integers(0, min(13, dimension + 1)))
            left_support = tuple(
                sorted(
                    int(value)
                    for value in rng.choice(dimension, left_size, replace=False)
                )
            )
            right_support = tuple(
                sorted(
                    int(value)
                    for value in rng.choice(dimension, right_size, replace=False)
                )
            )
            observed = moment(order, left_support, right_support)
            expected = enumerated_moment(
                left_values,
                right_values,
                left_support,
                right_support,
            )
            if observed != expected:
                raise AssertionError(
                    (order, left_support, right_support, observed, expected)
                )
            comparisons += 1

    # A vertical odd support and a horizontal odd support force exactly one
    # hidden-permutation match, hence magnitude 1/q at any order.
    vertical = (0, 32, 64, 96, 128)
    horizontal = (0, 1, 2)
    q32 = moment(32, vertical, horizontal)
    if abs(q32) != Fraction(1, 32):
        raise AssertionError(("q32 vertical-horizontal", q32))

    print(
        "signed-permutation link moment passed: "
        f"enumeration_comparisons={comparisons},q32_check={q32}"
    )


if __name__ == "__main__":
    main()
