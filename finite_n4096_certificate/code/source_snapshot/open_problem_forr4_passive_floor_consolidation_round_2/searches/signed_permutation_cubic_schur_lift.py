#!/usr/bin/env python3
"""Occurrence-split Schur lifts of the signed-permutation cubic endpoint."""

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


def balanced_split_formula(order: int) -> float:
    return (
        4
        * (order - 1)
        * (order * order - order + 1)
        / (order**4 * np.sqrt(order * order - 1))
    )


def normalized_nuclear(matrix: np.ndarray) -> tuple[float, np.ndarray]:
    rows, columns = matrix.shape
    if rows <= columns:
        gram = matrix @ matrix.T
    else:
        gram = matrix.T @ matrix
    eigenvalues = np.linalg.eigvalsh((gram + gram.T) / 2)
    eigenvalues = np.maximum(eigenvalues, 0)
    singular = np.sqrt(eigenvalues)
    return float(singular.sum() / np.sqrt(rows * columns)), singular


def endpoint_link_lift(
    moment: np.ndarray,
    degree_three: list[tuple[int, ...]],
    ket_size: int,
) -> tuple[np.ndarray, int, int]:
    bra_size = 3 - ket_size
    ket_supports = list(combinations(range(N), ket_size))
    bra_supports = list(combinations(range(N), bra_size))
    support_index = {support: index for index, support in enumerate(degree_three)}
    lift = np.zeros((len(ket_supports), len(bra_supports), N))
    for ket_index, ket in enumerate(ket_supports):
        ket_set = set(ket)
        for bra_index, bra in enumerate(bra_supports):
            if ket_set.intersection(bra):
                continue
            union = tuple(sorted(ket + bra))
            lift[ket_index, bra_index] = moment[support_index[union]]
    return lift, len(ket_supports), len(bra_supports)


def main() -> None:
    left_values, right_values = all_pairs()
    degree_one = supports(1)
    degree_three = supports(3)
    left_one = features(left_values, degree_one)
    left_three = features(left_values, degree_three)
    right_one = features(right_values, degree_one)
    m31 = left_three.T.astype(float) @ right_one / len(left_values)
    m11 = left_one.T.astype(float) @ right_one / len(left_values)

    rows = []
    worst = (0.0, None)
    worst_singular = np.array([])
    for ket_size in range(4):
        lift, ket_count, bra_count = endpoint_link_lift(
            m31, degree_three, ket_size
        )
        tensor = (
            lift[:, :, :, None, None]
            * m11[None, None, :, :, None]
            * m11[None, None, None, :, :]
        )
        # mask chooses which of singleton blocks 2,3,4 join the ket side.
        for mask in range(8):
            row_axes = [0] + [
                2 + block for block in range(3) if (mask >> block) & 1
            ]
            column_axes = [1] + [
                2 + block
                for block in range(3)
                if not ((mask >> block) & 1)
            ]
            matrix = np.transpose(tensor, row_axes + column_axes).reshape(
                ket_count * N ** mask.bit_count(),
                bra_count * N ** (3 - mask.bit_count()),
            )
            value, singular = normalized_nuclear(matrix)
            rows.append(
                f"endpoint_split={ket_size}|{3-ket_size},"
                f"singleton_mask={mask:03b},shape={matrix.shape[0]}x"
                f"{matrix.shape[1]},normalized_nuclear={value:.12g}"
            )
            if value > worst[0]:
                worst = (value, (ket_size, mask, matrix.shape))
                worst_singular = singular

    if worst[0] > 1 / Q + 3e-11:
        raise AssertionError(("cubic Schur lift exceeds minimal scale", worst))
    expected_worst = balanced_split_formula(Q)
    if not np.isclose(worst[0], expected_worst, atol=3e-10):
        raise AssertionError(("balanced split formula", worst, expected_worst))
    positive = worst_singular[worst_singular > 2e-10]
    rounded, counts = np.unique(np.round(positive, 10), return_counts=True)
    spectrum = ",".join(
        f"{value:.10g}x{count}" for value, count in zip(rounded, counts)
    )
    print(
        "signed-permutation cubic Schur lift:\n"
        + "\n".join(rows)
        + "\n"
        + f"worst={worst[0]:.12g},placement={worst[1]},minimal_scale={1/Q:.12g},"
        + f"rank={len(positive)},spectrum={spectrum},"
        + f"q32_formula={balanced_split_formula(32):.12g}"
    )


if __name__ == "__main__":
    main()
