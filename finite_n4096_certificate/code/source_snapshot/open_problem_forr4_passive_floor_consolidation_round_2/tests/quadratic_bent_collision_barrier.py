#!/usr/bin/env python3
"""Exact affine-flat collision barrier for quadratic-bent features."""

from __future__ import annotations

import random
from fractions import Fraction


SEED = 2026071419


def plane() -> frozenset[int]:
    return frozenset((0, 1, 2, 3))


def extensions(dimension: int) -> list[frozenset[int]]:
    base = plane()
    # One representative for each nonzero vector of V/base.  With the
    # chosen coordinate plane, representatives are multiples of four.
    result = []
    for representative in range(4, 1 << dimension, 4):
        flat = frozenset(x ^ shift for x in base for shift in (0, representative))
        complement = flat - base
        if len(flat) != 8 or len(complement) != 4:
            raise AssertionError(("affine 3-flat", dimension, representative, flat))
        result.append(frozenset(complement))
    if len(set(result)) != (1 << (dimension - 2)) - 1:
        raise AssertionError(("extension count", dimension, len(set(result))))
    return result


def quadratic_value(
    point: int,
    dimension: int,
    linear: int,
    constant: int,
    quadratic_mask: int,
) -> int:
    value = constant ^ ((point & linear).bit_count() & 1)
    edge = 0
    for first in range(dimension):
        for second in range(first + 1, dimension):
            if (quadratic_mask >> edge) & 1:
                value ^= ((point >> first) & 1) & ((point >> second) & 1)
            edge += 1
    return value


def feature(
    support: frozenset[int],
    dimension: int,
    linear: int,
    constant: int,
    quadratic_mask: int,
) -> int:
    parity = 0
    for point in support:
        parity ^= quadratic_value(
            point, dimension, linear, constant, quadratic_mask
        )
    return 1 if parity == 0 else -1


def invert_binary_matrix(rows: list[int], dimension: int) -> list[int] | None:
    augmented = [
        rows[row] | (1 << (dimension + row)) for row in range(dimension)
    ]
    for column in range(dimension):
        pivot = next(
            (
                row
                for row in range(column, dimension)
                if (augmented[row] >> column) & 1
            ),
            None,
        )
        if pivot is None:
            return None
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        for row in range(dimension):
            if row != column and ((augmented[row] >> column) & 1):
                augmented[row] ^= augmented[column]
    mask = (1 << dimension) - 1
    return [
        (augmented[row] >> dimension) & mask for row in range(dimension)
    ]


def canonical_plane_cross_moment(dimension: int) -> Fraction:
    edges = [
        (first, second)
        for first in range(dimension)
        for second in range(first + 1, dimension)
    ]
    total = 0
    signed = 0
    for quadratic_mask in range(1 << len(edges)):
        rows = [0] * dimension
        for edge, (first, second) in enumerate(edges):
            if (quadratic_mask >> edge) & 1:
                rows[first] |= 1 << second
                rows[second] |= 1 << first
        inverse = invert_binary_matrix(rows, dimension)
        if inverse is None:
            continue
        total += 1
        primal = (rows[0] >> 1) & 1
        dual = (inverse[0] >> 1) & 1
        signed += 1 if primal == dual else -1
    return Fraction(signed, total)


def gaussian_binomial(dimension: int, rank: int) -> int:
    result = Fraction(1)
    for offset in range(rank):
        result *= Fraction(
            (1 << (dimension - offset)) - 1,
            (1 << (rank - offset)) - 1,
        )
    if result.denominator != 1:
        raise AssertionError(("Gaussian binomial", dimension, rank, result))
    return result.numerator


def symplectic_fourier(dimension: int, half_rank: int) -> Fraction:
    result = Fraction((-1) ** half_rank)
    for step in range(1, half_rank + 1):
        result /= (1 << (dimension - 2 * step + 1)) - 1
    return result


def endpoint_five_one_norm_squared(dimension: int) -> Fraction:
    # For a five-set A with fixed XOR, S=A triangle {XOR(A)} is either:
    # - four nonzero zero-XOR points spanning dimension three (7 orbits
    #   per 3-space), whose alternating character has rank two; or
    # - five nonzero zero-XOR points spanning dimension four (168 orbits
    #   per 4-space), whose alternating character has rank four.
    rank_two = symplectic_fourier(dimension, 1)
    rank_four = symplectic_fourier(dimension, 2)
    return (
        7 * gaussian_binomial(dimension, 3) * rank_two**2
        + 168 * gaussian_binomial(dimension, 4) * rank_four**2
    )


def main() -> None:
    rng = random.Random(SEED)
    rows = []
    for dimension in (4, 6, 10):
        base = plane()
        complements = extensions(dimension)
        for _ in range(200):
            linear = rng.randrange(1 << dimension)
            constant = rng.randrange(2)
            quadratic_mask = rng.randrange(1 << (dimension * (dimension - 1) // 2))
            target = feature(
                base, dimension, linear, constant, quadratic_mask
            )
            if any(
                feature(
                    support,
                    dimension,
                    linear,
                    constant,
                    quadratic_mask,
                )
                != target
                for support in complements
            ):
                raise AssertionError(("quadratic feature collision", dimension))
        class_size = 1 + len(complements)
        expected = 1 << (dimension - 2)
        if class_size != expected:
            raise AssertionError(("collision class size", class_size, expected))
        rows.append(
            f"N={1 << dimension}:identical_degree4_columns={class_size}"
        )

    cross_rows = []
    for dimension in (4, 6):
        observed = canonical_plane_cross_moment(dimension)
        order = 1 << dimension
        expected = Fraction(2, order - 2)
        if observed != expected:
            raise AssertionError(
                ("canonical plane cross moment", dimension, observed, expected)
            )
        amplified = Fraction(order, 4) * observed
        cross_rows.append(
            f"N={order}:canonical_class_cross={observed},"
            f"class_amplification={amplified}"
        )

    endpoint_rows = []
    expected_n16 = Fraction(39, 7)
    if endpoint_five_one_norm_squared(4) != expected_n16:
        raise AssertionError(
            ("N=16 M51 formula", endpoint_five_one_norm_squared(4))
        )
    for dimension in (4, 6, 10):
        value = endpoint_five_one_norm_squared(dimension)
        endpoint_rows.append(
            f"N={1 << dimension}:M51_squared={value},"
            f"M51={float(value) ** 0.5:.12g}"
        )
    n1024 = endpoint_five_one_norm_squared(10)
    if n1024 != Fraction(11_182_413, 64_897):
        raise AssertionError(("N=1024 M51 value", n1024))

    print(
        "quadratic-bent collision barrier passed: "
        + ", ".join(rows + cross_rows + endpoint_rows)
    )


if __name__ == "__main__":
    main()
