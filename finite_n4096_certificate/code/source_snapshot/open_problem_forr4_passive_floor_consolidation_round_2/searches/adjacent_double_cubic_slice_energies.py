#!/usr/bin/env python3
"""Exact q=4 unordered split energies for the (3,3,1,1) record-one sector."""

from __future__ import annotations

from itertools import combinations

import numpy as np

from signed_permutation_full_sector_spectra import (
    N,
    Q,
    all_pairs,
    features,
    record_size,
    supports,
)


def incidence(
    parts: list[tuple[int, ...]],
    full_supports: list[tuple[int, ...]],
) -> np.ndarray:
    result = np.zeros((len(parts), len(full_supports)))
    part_size = len(parts[0]) if parts else 0
    if part_size == 0:
        result[0] = 1
        return result
    part_index = {part: index for index, part in enumerate(parts)}
    for column, support in enumerate(full_supports):
        for part in combinations(support, part_size):
            result[part_index[part], column] = 1
    return result


def main() -> None:
    left_values, right_values = all_pairs()
    degree_one = supports(1)
    degree_three = supports(3)
    left_three = features(left_values, degree_three)
    right_one = features(right_values, degree_one)
    right_three = features(right_values, degree_three)
    normalization = len(left_values)
    m33 = left_three.T.astype(float) @ right_three / normalization
    m31 = left_three.T.astype(float) @ right_one / normalization

    row_record = np.array(
        [record_size(support, "right") for support in degree_three]
    )
    column_record = np.array(
        [record_size(support, "left") for support in degree_three]
    )
    endpoint = column_record == 1
    l_shape = (row_record == 1) & (column_record == 1)
    endpoint_supports = [
        support for support, keep in zip(degree_three, endpoint) if keep
    ]
    l_supports = [
        support for support, keep in zip(degree_three, l_shape) if keep
    ]
    squared = (
        np.abs(m33[np.ix_(endpoint, l_shape)])[:, :, None] ** 2
        * np.abs(m31[l_shape, :])[None, :, :] ** 2
    )
    summed_singleton = squared.sum(axis=2)

    endpoint_incidence = {}
    l_incidence = {}
    for size in range(4):
        endpoint_parts = list(combinations(range(N), size))
        l_parts = list(combinations(range(N), size))
        endpoint_incidence[size] = incidence(endpoint_parts, endpoint_supports)
        l_incidence[size] = incidence(l_parts, l_supports)

    rows = []
    worst = (0.0, None)
    for endpoint_ket in range(4):
        for middle_ket in range(4):
            row_energy_matrix = (
                endpoint_incidence[endpoint_ket]
                @ summed_singleton
                @ l_incidence[middle_ket].T
            )
            row_energy = float(np.max(row_energy_matrix))

            endpoint_bra = 3 - endpoint_ket
            middle_bra = 3 - middle_ket
            column_energy = 0.0
            left_incidence = endpoint_incidence[endpoint_bra]
            right_incidence = l_incidence[middle_bra]
            for singleton in range(N):
                value = (
                    left_incidence
                    @ squared[:, :, singleton]
                    @ right_incidence.T
                )
                column_energy = max(column_energy, float(np.max(value)))

            coefficient = min(np.sqrt(row_energy), np.sqrt(column_energy))
            rows.append(
                f"endpoint_split={endpoint_ket}|{3-endpoint_ket},"
                f"middle_split={middle_ket}|{3-middle_ket},"
                f"row_energy={row_energy:.12g},"
                f"column_energy={column_energy:.12g},"
                f"coefficient={coefficient:.12g}"
            )
            if coefficient > worst[0]:
                worst = (
                    coefficient,
                    (endpoint_ket, middle_ket, row_energy, column_energy),
                )

    print(
        "adjacent double-cubic slice energies:\n"
        + "\n".join(rows)
        + "\n"
        + f"worst={worst[0]:.12g},placement={worst[1]}"
    )


if __name__ == "__main__":
    main()
