#!/usr/bin/env python3
"""Exact entry bounds for the three easier double-cubic profiles."""

from __future__ import annotations

from math import comb

import numpy as np

from signed_permutation_full_sector_spectra import (
    Q,
    all_pairs,
    features,
    record_size,
    supports,
)


def main() -> None:
    left_values, right_values = all_pairs()
    degree_one = supports(1)
    degree_three = supports(3)
    left_one = features(left_values, degree_one)
    left_three = features(left_values, degree_three)
    right_one = features(right_values, degree_one)
    right_three = features(right_values, degree_three)
    normalization = len(left_values)
    m13 = left_one.T.astype(float) @ right_three / normalization
    m31 = left_three.T.astype(float) @ right_one / normalization
    m33 = left_three.T.astype(float) @ right_three / normalization

    row_record = np.array(
        [record_size(support, "right") for support in degree_three]
    )
    column_record = np.array(
        [record_size(support, "left") for support in degree_three]
    )

    # Profile (3,1,3,1): the second cubic block is compatible with
    # singleton records on both sides and is therefore an L-shape.
    l_shape = (row_record == 1) & (column_record == 1)
    endpoint_record = column_record == 1
    endpoint_link = np.max(np.abs(m31[endpoint_record]))
    middle_left = np.max(np.abs(m13[:, l_shape]))
    middle_right = np.max(np.abs(m31[l_shape, :]))
    separated_product = float(endpoint_link * middle_left * middle_right)
    expected_record_one = 1 / (Q**3 * (Q - 1) ** 2)
    if not np.isclose(separated_product, expected_record_one, atol=2e-12):
        raise AssertionError(
            ("separated double-cubic entry", separated_product)
        )

    # Profile (1,3,3,1), central record one: both cubic supports are L-shapes.
    central_one = m33[np.ix_(l_shape, l_shape)]
    adjacent_one = float(
        np.max(np.abs(m13[:, l_shape]))
        * np.max(np.abs(central_one))
        * np.max(np.abs(m31[l_shape, :]))
    )
    if not np.isclose(adjacent_one, expected_record_one, atol=2e-12):
        raise AssertionError(("adjacent cubic record one", adjacent_one))

    # Central record three: the outer singleton links can have full 1/q
    # coherence, while the central match costs 1/C(q,3).
    left_three_record = (row_record == 1) & (column_record == 3)
    right_three_record = (row_record == 3) & (column_record == 1)
    central_three = m33[np.ix_(left_three_record, right_three_record)]
    adjacent_three = float(
        np.max(np.abs(m13[:, left_three_record]))
        * np.max(np.abs(central_three))
        * np.max(np.abs(m31[right_three_record, :]))
    )
    expected_record_three = 1 / (Q**2 * comb(Q, 3))
    if not np.isclose(adjacent_three, expected_record_three, atol=2e-12):
        raise AssertionError(("adjacent cubic record three", adjacent_three))

    q32_record_one = 32 / 31**2
    q32_record_three = 32**2 / comb(32, 3)
    endpoint_to_l = m33[np.ix_(endpoint_record, l_shape)]
    endpoint_l_max_row_energy = float(
        np.max(np.sum(np.abs(endpoint_to_l) ** 2, axis=1))
    )
    endpoint_l_max_column_energy = float(
        np.max(np.sum(np.abs(endpoint_to_l) ** 2, axis=0))
    )
    endpoint_l_op = float(
        np.linalg.svd(endpoint_to_l, compute_uv=False)[0]
    )
    expected_endpoint_l_energy = 1 + 2 * Q / (Q - 1) ** 2
    if not np.isclose(
        endpoint_l_max_row_energy,
        expected_endpoint_l_energy,
        atol=2e-12,
    ):
        raise AssertionError(
            ("endpoint-to-L row energy", endpoint_l_max_row_energy)
        )
    q32_endpoint_l_alternating = (
        np.sqrt(1 + 64 / 31**2) / (32 * 31)
    )
    print(
        "signed-permutation double-cubic entries passed: "
        f"q={Q},L_shapes={int(l_shape.sum())},"
        f"separated_record_one={separated_product:.12g},"
        f"adjacent_record_one={adjacent_one:.12g},"
        f"adjacent_record_three={adjacent_three:.12g},"
        f"endpoint_L_max_row_energy={endpoint_l_max_row_energy:.12g},"
        "endpoint_L_max_column_energy="
        f"{endpoint_l_max_column_energy:.12g},"
        f"endpoint_L_op={endpoint_l_op:.12g},"
        f"q32_rank_record_one={q32_record_one:.12g},"
        f"q32_rank_record_three={q32_record_three:.12g},"
        "q32_endpoint_L_alternating="
        f"{q32_endpoint_l_alternating:.12g}"
    )


if __name__ == "__main__":
    main()
