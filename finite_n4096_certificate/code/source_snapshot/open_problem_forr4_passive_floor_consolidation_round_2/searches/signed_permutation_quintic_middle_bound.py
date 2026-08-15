#!/usr/bin/env python3
"""Exact q=4 check and q=32 bound for a quintic middle decoration."""

from __future__ import annotations

import numpy as np

from signed_permutation_full_sector_spectra import (
    N,
    Q,
    all_pairs,
    features,
    record_size,
    supports,
)


def middle_coefficient(order: int) -> float:
    return order / (order - 1) ** 2


def main() -> None:
    left_values, right_values = all_pairs()
    degree_one = supports(1)
    degree_five = supports(5)
    left_one = features(left_values, degree_one)
    left_five = features(left_values, degree_five)
    right_one = features(right_values, degree_one)
    right_five = features(right_values, degree_five)
    m15 = left_one.T.astype(float) @ right_five / len(left_values)
    m51 = left_five.T.astype(float) @ right_one / len(left_values)

    compatible = np.array(
        [
            record_size(support, "left") == 1
            and record_size(support, "right") == 1
            for support in degree_five
        ]
    )
    compatible_count = int(compatible.sum())
    restricted15 = m15[:, compatible]
    restricted51 = m51[compatible, :]
    maximum15 = float(np.max(np.abs(restricted15)))
    maximum51 = float(np.max(np.abs(restricted51)))
    expected_link = 1 / (Q * (Q - 1))
    if not np.isclose(maximum15, expected_link, atol=2e-12):
        raise AssertionError(("M15 compatible coherence", maximum15))
    if not np.isclose(maximum51, expected_link, atol=2e-12):
        raise AssertionError(("M51 compatible coherence", maximum51))

    adjacent_product = np.abs(restricted15[:, :, None]) * np.abs(
        restricted51[None, :, :]
    )
    maximum_product = float(np.max(adjacent_product))
    expected_product = expected_link**2
    if not np.isclose(maximum_product, expected_product, atol=2e-12):
        raise AssertionError(
            ("quintic adjacent product coherence", maximum_product)
        )

    q32 = middle_coefficient(32)
    if q32 >= 1 / 30:
        raise AssertionError(("q=32 quintic coefficient", q32))
    print(
        "signed-permutation quintic middle bound passed: "
        f"q={Q},supports={len(degree_five)},compatible={compatible_count},"
        f"max_M15={maximum15:.12g},max_M51={maximum51:.12g},"
        f"max_adjacent_product={maximum_product:.12g},"
        f"q32_rank_Frobenius={q32:.12g},"
        f"q32_relative_to_1_over_q={32*q32:.12g}"
    )


if __name__ == "__main__":
    main()
