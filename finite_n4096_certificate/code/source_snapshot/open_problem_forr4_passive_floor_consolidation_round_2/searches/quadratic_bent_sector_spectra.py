#!/usr/bin/env python3
"""Low-copy link spectra for the complete N=16 quadratic-bent orbit."""

from __future__ import annotations

from itertools import combinations

import numpy as np


N = 16
MAX_DEGREE = 6


def sylvester(order: int) -> np.ndarray:
    result = np.array([[1]], dtype=int)
    while result.shape[0] < order:
        result = np.block([[result, result], [result, -result]])
    return result


def all_bent_pairs() -> tuple[np.ndarray, np.ndarray]:
    hadamard = sylvester(N)
    x_rows = []
    y_rows = []
    coordinates = np.arange(N)
    for mask in range(1 << N):
        x = (1 - 2 * ((mask >> coordinates) & 1)).astype(np.int8)
        transform = hadamard @ x
        if np.all(np.abs(transform) == np.sqrt(N)):
            x_rows.append(x)
            y_rows.append((transform // int(np.sqrt(N))).astype(np.int8))
    return np.array(x_rows), np.array(y_rows)


def feature_matrix(values: np.ndarray, degree: int) -> np.ndarray:
    supports = list(combinations(range(N), degree))
    result = np.empty((len(values), len(supports)), dtype=np.int8)
    for column, support in enumerate(supports):
        result[:, column] = np.prod(values[:, support], axis=1)
    return result


def low_rank_singular_norm(
    left: np.ndarray, right: np.ndarray
) -> float:
    count = len(left)
    left_gram = left.astype(float) @ left.T / count
    right_gram = right.astype(float) @ right.T / count
    left_values, left_vectors = np.linalg.eigh(left_gram)
    right_values, right_vectors = np.linalg.eigh(right_gram)
    left_keep = left_values > 1e-9
    right_keep = right_values > 1e-9
    core = (
        np.sqrt(left_values[left_keep])[:, None]
        * (left_vectors[:, left_keep].T @ right_vectors[:, right_keep])
        * np.sqrt(right_values[right_keep])[None, :]
    )
    return float(np.linalg.svd(core, compute_uv=False)[0])


def main() -> None:
    x_values, y_values = all_bent_pairs()
    if len(x_values) != 896:
        raise AssertionError(("N=16 bent count", len(x_values)))
    x_features = [
        feature_matrix(x_values, degree)
        for degree in range(1, MAX_DEGREE + 1)
    ]
    y_features = [
        feature_matrix(y_values, degree)
        for degree in range(1, MAX_DEGREE + 1)
    ]

    # Parity forces the checkerboard zeros.  The nonzero values are the
    # independently recorded exact N=16 spectrum, rounded only at 1e-6.
    expected = np.array(
        [
            [1, 0, 0.845154254729, 0, 2.360387377409, 0],
            [0, 1.142857142857, 0, 1.616244071284, 0, 4.426266681731],
            [0.845154254729, 0, 1.714285714286, 0, 3.833259389999, 0],
            [0, 1.616244071284, 0, 6.857142857143, 0, 9.142857142857],
            [2.360387377409, 0, 3.833259389999, 0, 13.714285714286, 0],
            [0, 4.426266681731, 0, 9.142857142857, 0, 36.571428571429],
        ]
    )
    observed = np.empty_like(expected)
    for left in range(MAX_DEGREE):
        for right in range(MAX_DEGREE):
            observed[left, right] = low_rank_singular_norm(
                x_features[left], y_features[right]
            )
    if not np.allclose(observed, expected, atol=2e-6):
        raise AssertionError(("quadratic-bent sector spectrum", observed))

    signed_permutation_diagonal = np.array(
        [1, 8 / 3, 4, 16, 32, 256 / 3]
    )
    if not np.all(np.diag(observed) < signed_permutation_diagonal + 1e-12):
        raise AssertionError(("diagonal comparison", np.diag(observed)))

    print(
        "quadratic-bent sector spectra passed: "
        f"bent_functions={len(x_values)}, "
        "diagonal="
        + ",".join(f"{value:.12g}" for value in np.diag(observed))
        + ", signed_permutation_diagonal="
        + ",".join(f"{value:.12g}" for value in signed_permutation_diagonal)
    )


if __name__ == "__main__":
    main()
