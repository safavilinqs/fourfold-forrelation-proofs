#!/usr/bin/env python3
"""Exact counterexample to the proposed near-unit temporal square function."""

from __future__ import annotations

from fractions import Fraction
from itertools import product


K4 = (
    (1, 1, 1, 1),
    (1, -1, 1, -1),
    (1, 1, -1, -1),
    (1, -1, -1, 1),
)


def probability(outcome: int, signs: tuple[int, ...]) -> Fraction:
    amplitude_numerator = sum(
        K4[outcome][block] * signs[block] for block in range(4)
    )
    return Fraction(amplitude_numerator**2, 16)


def signed_output_matrix() -> tuple[tuple[Fraction, ...], ...]:
    result = []
    for first in range(4):
        row = []
        for second in range(4):
            value = Fraction()
            for signs in product((-1, 1), repeat=4):
                character = signs[0] * signs[1] * signs[2] * signs[3]
                value += (
                    Fraction(character, 16)
                    * probability(first, signs)
                    * probability(second, signs)
                )
            row.append(value)
        result.append(tuple(row))
    return tuple(result)


def main() -> None:
    matrix = signed_output_matrix()
    expected = tuple(
        tuple(Fraction(3 if row == column else -1, 32) for column in range(4))
        for row in range(4)
    )
    if matrix != expected:
        raise AssertionError(("signed output matrix", matrix, expected))

    transcript_mass = sum(abs(value) for row in matrix for value in row)
    if transcript_mass != Fraction(3, 4):
        raise AssertionError(("transcript mass", transcript_mass))

    # A dose-one probe is uniform on the four one-photon block supports.
    # For each two-block mask, the local ket/bra square mass is 1/8.
    # With two nodes there are six ordered complementary masks.
    joint_square_mass = 6 * Fraction(1, 8) ** 2
    ratio_squared_n1 = transcript_mass**2 / joint_square_mass
    if joint_square_mass != Fraction(3, 32) or ratio_squared_n1 != 6:
        raise AssertionError(
            ("joint square comparison", joint_square_mass, ratio_squared_n1)
        )

    # Restrict an N=2 chain to coordinate zero in every block.  Its graph
    # entry is 2^{-3/2}, whereas the proposed common graph factor is
    # 2^{-1/2}.  The squared required frame constant is therefore 6/2^2.
    ratio_squared_n2 = ratio_squared_n1 / 4
    threshold_squared = Fraction(32, 27) ** 2
    if ratio_squared_n2 <= threshold_squared:
        raise AssertionError(
            ("counterexample no longer crosses target", ratio_squared_n2)
        )

    print(
        "two-copy square-function target falsified: "
        f"transcript_mass={transcript_mass}, joint_square_mass={joint_square_mass}, "
        f"N1_ratio_squared={ratio_squared_n1}, "
        f"N2_ratio_squared={ratio_squared_n2}, "
        f"target_squared={threshold_squared}"
    )


if __name__ == "__main__":
    main()
