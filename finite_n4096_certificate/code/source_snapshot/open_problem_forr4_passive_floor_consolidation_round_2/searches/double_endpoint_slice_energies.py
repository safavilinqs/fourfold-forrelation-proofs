#!/usr/bin/env python3
"""Exact cubic endpoint slice energies and the two-endpoint lift bound."""

from __future__ import annotations

from itertools import combinations

import numpy as np

from signed_permutation_full_sector_spectra import (
    N,
    Q,
    all_pairs,
    features,
    supports,
)


def slice_energies(order: int) -> tuple[float, float, float, float]:
    return (
        (order * order + 2) / 6,
        (order * order + 2) / (2 * order * order),
        (order * order - 2 * order + 2)
        / (order * order * (order - 1)),
        1 / (order * order),
    )


def two_endpoint_coefficient(order: int) -> float:
    _, energy_one, energy_two, _ = slice_energies(order)
    return np.sqrt(energy_one * energy_two)


def main() -> None:
    left_values, right_values = all_pairs()
    degree_one = supports(1)
    degree_three = supports(3)
    left_three = features(left_values, degree_three)
    right_one = features(right_values, degree_one)
    m31 = left_three.T.astype(float) @ right_one / len(left_values)
    support_index = {support: index for index, support in enumerate(degree_three)}

    observed = []
    # k=0: the neighboring singleton is fixed and all cubic marks are summed.
    observed.append(float(np.max(np.sum(np.abs(m31) ** 2, axis=0))))

    for ket_size in (1, 2):
        fixed_parts = list(combinations(range(N), ket_size))
        summed_parts = list(combinations(range(N), 3 - ket_size))
        maximum = 0.0
        for fixed in fixed_parts:
            fixed_set = set(fixed)
            for singleton in range(N):
                energy = 0.0
                for summed in summed_parts:
                    if fixed_set.intersection(summed):
                        continue
                    union = tuple(sorted(fixed + summed))
                    energy += abs(m31[support_index[union], singleton]) ** 2
                maximum = max(maximum, energy)
        observed.append(maximum)

    observed.append(float(np.max(np.abs(m31)) ** 2))
    expected = slice_energies(Q)
    if not np.allclose(observed, expected, atol=2e-12):
        raise AssertionError(("cubic endpoint slice energies", observed, expected))

    # Exhaust the four-by-four split table algebraically.  Whole-block
    # complementary placements use the sharper block-coherent 1/q bound;
    # the internal maximum must occur at a 1|2 pairing.
    worst = 0.0
    placement = None
    for left_split in range(4):
        for right_split in range(4):
            row = np.sqrt(expected[left_split] * expected[right_split])
            column = np.sqrt(
                expected[3 - left_split] * expected[3 - right_split]
            )
            value = min(row, column)
            if {left_split, right_split} == {0, 3}:
                value = min(value, 1 / Q)
            if value > worst:
                worst = value
                placement = (left_split, right_split)
    target = two_endpoint_coefficient(Q)
    if not np.isclose(worst, target, atol=2e-12):
        raise AssertionError(("two-endpoint split maximum", worst, target))

    q32 = two_endpoint_coefficient(32)
    beta = 5 / 6
    print(
        "double-endpoint slice energies passed: "
        f"q={Q},E0={observed[0]:.12g},E1={observed[1]:.12g},"
        f"E2={observed[2]:.12g},E3={observed[3]:.12g},"
        f"worst_split={placement},q4_coefficient={worst:.12g},"
        f"q32_coefficient={q32:.12g},"
        f"q32_beta8_coefficient={q32*beta**8:.12g}"
    )


if __name__ == "__main__":
    main()
