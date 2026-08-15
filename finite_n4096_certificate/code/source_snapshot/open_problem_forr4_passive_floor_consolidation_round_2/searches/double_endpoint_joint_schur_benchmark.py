#!/usr/bin/env python3
"""Joint Schur benchmarks for the two cubic endpoint profile."""

from __future__ import annotations

from itertools import combinations, permutations, product

import numpy as np


def sylvester(order: int) -> np.ndarray:
    result = np.array([[1]], dtype=np.int8)
    while result.shape[0] < order:
        result = np.block([[result, result], [result, -result]])
    return result


def planted_values(order: int) -> tuple[np.ndarray, np.ndarray]:
    walsh = sylvester(order)
    left = []
    right = []
    for permutation in permutations(range(order)):
        for signs in product((-1, 1), repeat=order):
            signed_permutation = np.zeros((order, order), dtype=np.int8)
            for column, row in enumerate(permutation):
                signed_permutation[row, column] = signs[column]
            left.append((walsh @ signed_permutation).reshape(-1))
            right.append((signed_permutation @ walsh).reshape(-1))
    return np.array(left), np.array(right)


def endpoint_slice(order: int) -> tuple[np.ndarray, int]:
    dimension = order * order
    left, right = planted_values(order)
    triples = list(combinations(range(dimension), 3))
    pairs = list(combinations(range(dimension), 2))
    left_features = np.array(
        [[np.prod(value[list(support)]) for support in triples] for value in left]
    )
    moment = left_features.T.astype(float) @ right / len(left)
    triple_index = {support: index for index, support in enumerate(triples)}
    tensor = np.zeros((dimension, dimension, len(pairs)))
    for singleton in range(dimension):
        for pair_index, pair in enumerate(pairs):
            if singleton in pair:
                continue
            support = tuple(sorted((singleton,) + pair))
            tensor[singleton, :, pair_index] = moment[triple_index[support]]
    return tensor, len(pairs)


def same_side_value(order: int, endpoint: np.ndarray, pair_count: int):
    dimension = order * order
    frame = endpoint.reshape(dimension, dimension * pair_count)
    gram = frame @ frame.T
    expected_energy = (order * order + 2) / 2
    if not np.allclose(
        gram, expected_energy * np.eye(dimension), atol=3e-11
    ):
        raise AssertionError(("endpoint tight frame", order, gram))
    value = expected_energy / (order * (dimension * (dimension - 1) / 2))
    expected_formula = (order * order + 2) / (
        order**3 * (order * order - 1)
    )
    if not np.isclose(value, expected_formula, atol=2e-14):
        raise AssertionError(("same-side formula", order, value))
    return value, expected_energy


def alternating_value(order: int, endpoint: np.ndarray, pair_count: int):
    dimension = order * order
    walsh = sylvester(order).astype(float)
    hadamard = np.kron(walsh, walsh) / order

    # Rows are (left endpoint singleton, first middle singleton,
    # right endpoint singleton); columns contain the two endpoint pairs
    # and the second middle singleton.  Form the row Gram in N-by-N
    # middle-label blocks, avoiding the much wider rectangular matrix.
    left_gram = np.einsum(
        "ibe,jfe->bfij", endpoint, endpoint, optimize=True
    )
    right_gram = np.einsum(
        "bc,fc,kce,lce->bfkl",
        hadamard,
        hadamard,
        endpoint,
        endpoint,
        optimize=True,
    )
    row_count = dimension**3
    gram = np.empty((row_count, row_count))
    for first_middle in range(dimension):
        for right_endpoint in range(dimension):
            row_slice = slice(
                (first_middle * dimension + right_endpoint) * dimension,
                (first_middle * dimension + right_endpoint + 1) * dimension,
            )
            for other_middle in range(dimension):
                for other_endpoint in range(dimension):
                    column_slice = slice(
                        (other_middle * dimension + other_endpoint) * dimension,
                        (other_middle * dimension + other_endpoint + 1)
                        * dimension,
                    )
                    gram[row_slice, column_slice] = (
                        left_gram[first_middle, other_middle]
                        * right_gram[
                            first_middle,
                            other_middle,
                            right_endpoint,
                            other_endpoint,
                        ]
                    )
    eigenvalues = np.linalg.eigvalsh((gram + gram.T) / 2)
    singular_sum = float(np.sqrt(np.maximum(eigenvalues, 0)).sum())
    column_count = pair_count * pair_count * dimension
    return singular_sum / np.sqrt(row_count * column_count), int(
        np.sum(eigenvalues > 1e-10)
    )


def main() -> None:
    expected_alternating = {
        2: 0.47159181589114324,
        4: 0.06420087162467479,
    }
    rows = []
    for order in (2, 4):
        endpoint, pair_count = endpoint_slice(order)
        same_side, energy = same_side_value(order, endpoint, pair_count)
        alternating, rank = alternating_value(order, endpoint, pair_count)
        if not np.isclose(
            alternating, expected_alternating[order], atol=3e-11
        ):
            raise AssertionError(
                ("alternating joint Schur value", order, alternating)
            )
        rows.append(
            f"q={order},N={order * order},endpoint_frame_energy={energy:.12g},"
            f"same_side={same_side:.12g},alternating={alternating:.12g},"
            f"alternating_rank={rank}"
        )
    q32_same_side = (32 * 32 + 2) / (32**3 * (32 * 32 - 1))
    print(
        "double-endpoint joint Schur benchmark passed:\n"
        + "\n".join(rows)
        + "\n"
        + f"q32_same_side_formula={q32_same_side:.12g}"
    )


if __name__ == "__main__":
    main()
