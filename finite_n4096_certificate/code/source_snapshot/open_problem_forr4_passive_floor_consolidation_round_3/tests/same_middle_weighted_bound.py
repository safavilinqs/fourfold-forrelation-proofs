#!/usr/bin/env python3
"""Regression for the exact same-middle double-endpoint contraction."""

from __future__ import annotations

from itertools import combinations, permutations, product
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from alternating_double_endpoint_spectrum import scaled_endpoint_weights
from mixed_endpoint_weighted_bound import sylvester
from same_middle_weighted_bound import bound, deterministic_ledger


def endpoint_data(order: int, difference: int) -> tuple[float, float]:
    dimension = order * order
    left = np.array(
        [base for base in range(dimension) if base < (base ^ difference)],
        dtype=np.int32,
    )
    right = left ^ difference
    weights = np.stack(
        [
            scaled_endpoint_weights(order, singleton, left, right)
            for singleton in range(dimension)
        ]
    ).astype(float) / (order - 1)
    squared_column_energy = float(np.sum(weights[:, 0] ** 2))
    nuclear = float(np.linalg.svd(weights, compute_uv=False).sum())
    return squared_column_energy, nuclear


def direct_q2(
    joint_law: np.ndarray,
    row_difference_law: np.ndarray,
    column_difference_law: np.ndarray,
) -> tuple[float, float]:
    order = 2
    dimension = order * order
    pairs = list(combinations(range(dimension), 2))
    pair_left = np.array([pair[0] for pair in pairs], dtype=np.int32)
    pair_right = np.array([pair[1] for pair in pairs], dtype=np.int32)
    pair_difference = pair_left ^ pair_right
    pair_count = len(pairs)
    hadamard = sylvester(dimension) / np.sqrt(dimension)
    endpoint = np.empty((dimension, pair_count, dimension))
    for singleton in range(dimension):
        endpoint[singleton] = (
            scaled_endpoint_weights(
                order, singleton, pair_left, pair_right
            )[:, None]
            * hadamard[singleton ^ pair_difference]
        )
    tensor = np.einsum(
        "ieb,bc,dfc->iebcdf", endpoint, hadamard, endpoint
    )

    # Equal endpoint orientation: rows (i,d), columns (E,b,c,F).
    equal_kernel = np.transpose(tensor, (0, 4, 1, 2, 3, 5)).reshape(
        dimension**2, pair_count * dimension**2 * pair_count
    )
    equal_row_law = np.full(dimension**2, 1 / dimension**2)
    equal_column_law = np.array(
        [
            4
            * joint_law[
                int(pair_difference[left_pair]) - 1,
                int(pair_difference[right_pair]) - 1,
            ]
            / dimension**4
            for left_pair in range(pair_count)
            for first_middle in range(dimension)
            for second_middle in range(dimension)
            for right_pair in range(pair_count)
        ]
    )
    equal_value = float(
        np.linalg.svd(
            np.sqrt(equal_row_law)[:, None]
            * equal_kernel
            * np.sqrt(equal_column_law)[None, :],
            compute_uv=False,
        ).sum()
    )

    # Mixed endpoint orientation: rows (i,b,c,F), columns (E,d).
    mixed_kernel = np.transpose(tensor, (0, 2, 3, 5, 1, 4)).reshape(
        dimension**3 * pair_count, pair_count * dimension
    )
    mixed_row_law = np.array(
        [
            2
            * row_difference_law[int(pair_difference[pair]) - 1]
            / dimension**4
            for singleton in range(dimension)
            for first_middle in range(dimension)
            for second_middle in range(dimension)
            for pair in range(pair_count)
        ]
    )
    mixed_column_law = np.array(
        [
            2
            * column_difference_law[int(pair_difference[pair]) - 1]
            / dimension**2
            for pair in range(pair_count)
            for singleton in range(dimension)
        ]
    )
    mixed_value = float(
        np.linalg.svd(
            np.sqrt(mixed_row_law)[:, None]
            * mixed_kernel
            * np.sqrt(mixed_column_law)[None, :],
            compute_uv=False,
        ).sum()
    )
    return equal_value, mixed_value


def direct_hybrid_q2() -> tuple[float, float]:
    order = 2
    dimension = 4
    pairs = list(combinations(range(dimension), 2))
    pair_left = np.array([pair[0] for pair in pairs], dtype=np.int32)
    pair_right = np.array([pair[1] for pair in pairs], dtype=np.int32)
    pair_difference = pair_left ^ pair_right
    pair_count = len(pairs)
    hadamard = sylvester(dimension) / 2
    endpoint = np.empty((dimension, pair_count, dimension))
    for singleton in range(dimension):
        endpoint[singleton] = (
            scaled_endpoint_weights(
                order, singleton, pair_left, pair_right
            )[:, None]
            * hadamard[singleton ^ pair_difference]
        )

    small_walsh = sylvester(order)
    left_values = []
    right_values = []
    for permutation in permutations(range(order)):
        for signs in product((-1, 1), repeat=order):
            signed_permutation = np.zeros((order, order))
            for column, row in enumerate(permutation):
                signed_permutation[row, column] = signs[column]
            left_values.append((small_walsh @ signed_permutation).reshape(-1))
            right_values.append((signed_permutation @ small_walsh).reshape(-1))
    left_values = np.array(left_values)
    right_values = np.array(right_values)
    triples = list(combinations(range(dimension), 3))
    triple_features = np.array(
        [
            [np.prod(value[list(triple)]) for triple in triples]
            for value in right_values
        ]
    )
    whole_link = left_values.T @ triple_features / len(left_values)
    tensor = np.einsum(
        "ieb,bc,cd->iebcd", endpoint, hadamard, whole_link
    )

    difference_law = np.full(dimension - 1, 1 / (dimension - 1))
    pair_law = np.array(
        [
            2 * difference_law[int(pair_difference[pair]) - 1] / dimension
            for pair in range(pair_count)
        ]
    )

    # Rows (i,b,c,D), columns E.
    pair_column_kernel = np.transpose(tensor, (0, 2, 3, 4, 1)).reshape(
        dimension**3 * len(triples), pair_count
    )
    pair_column_row_law = np.full(
        pair_column_kernel.shape[0], 1 / pair_column_kernel.shape[0]
    )
    pair_column = float(
        np.linalg.svd(
            np.sqrt(pair_column_row_law)[:, None]
            * pair_column_kernel
            * np.sqrt(pair_law)[None, :],
            compute_uv=False,
        ).sum()
    )

    # Rows (i,D), columns (E,b,c).
    whole_column_kernel = np.transpose(tensor, (0, 4, 1, 2, 3)).reshape(
        dimension * len(triples), pair_count * dimension**2
    )
    whole_column_row_law = np.full(
        whole_column_kernel.shape[0], 1 / whole_column_kernel.shape[0]
    )
    whole_column_law = np.array(
        [
            pair_law[pair] / dimension**2
            for pair in range(pair_count)
            for first_middle in range(dimension)
            for second_middle in range(dimension)
        ]
    )
    whole_column = float(
        np.linalg.svd(
            np.sqrt(whole_column_row_law)[:, None]
            * whole_column_kernel
            * np.sqrt(whole_column_law)[None, :],
            compute_uv=False,
        ).sum()
    )
    return pair_column, whole_column


def main() -> None:
    for order in (4, 8, 16):
        large = order * order - 2 * order + 2
        vertical_energy, vertical_nuclear = endpoint_data(order, order)
        other_energy, other_nuclear = endpoint_data(order, 1)
        # This is the fixed-pair column energy.  It is twice the fixed-
        # singleton row energy because an XOR orbit contains N/2 pairs.
        if not np.isclose(
            vertical_energy, large / (order - 1), atol=2e-12
        ):
            raise AssertionError(("vertical endpoint energy", order))
        if not np.isclose(
            other_energy, 2 / (order - 1), atol=2e-12
        ):
            raise AssertionError(("other endpoint energy", order))
        if not np.isclose(
            vertical_nuclear, np.sqrt(2) * large, atol=2e-10
        ):
            raise AssertionError(("vertical endpoint nuclear", order))
        if not np.isclose(
            other_nuclear, np.sqrt(2) * order, atol=2e-10
        ):
            raise AssertionError(("other endpoint nuclear", order))

    rng = np.random.default_rng(2026071421)
    joint = rng.random((3, 3))
    joint /= joint.sum()
    row_law = rng.random(3)
    row_law /= row_law.sum()
    column_law = rng.random(3)
    column_law /= column_law.sum()
    direct_equal, direct_mixed = direct_q2(
        joint, row_law, column_law
    )
    energies = np.ones(3)
    endpoint_nuclear = np.full(3, 2 * np.sqrt(2))
    reduced_equal = 2 / 4**1.5 * np.sqrt(
        np.sum(joint * energies[:, None] * energies[None, :])
    )
    reduced_mixed = (
        2
        / 4**2.5
        * np.sqrt(row_law @ energies)
        * np.sum(np.sqrt(column_law) * endpoint_nuclear)
    )
    if not np.isclose(direct_equal, reduced_equal, atol=3e-13):
        raise AssertionError(("direct equal q2", direct_equal, reduced_equal))
    if not np.isclose(direct_mixed, reduced_mixed, atol=3e-13):
        raise AssertionError(("direct mixed q2", direct_mixed, reduced_mixed))

    hybrid_pair, hybrid_whole = direct_hybrid_q2()
    q2 = bound(2)
    if not np.isclose(
        hybrid_pair, float(q2.hybrid_pair_column), atol=3e-13
    ):
        raise AssertionError(("direct hybrid pair q2", hybrid_pair, q2))
    if not np.isclose(
        hybrid_whole, float(q2.hybrid_whole_column), atol=3e-13
    ):
        raise AssertionError(("direct hybrid whole q2", hybrid_whole, q2))

    q32 = bound(32)
    if not q32.equal_endpoint_orientation < q32.mixed_endpoint_orientation:
        raise AssertionError("q32 endpoint orientation ordering")
    raw, attenuated, margin = deterministic_ledger()
    slack = margin - attenuated
    if not 0.077 < float(slack) < 0.079:
        raise AssertionError(("double-endpoint finite-size slack", slack))
    print(
        "same-middle weighted contraction passed: "
        f"equal_q32={q32.equal_endpoint_orientation},"
        f"mixed_q32={q32.mixed_endpoint_orientation},"
        f"hybrid_q32={q32.hybrid_upper},"
        f"attenuated_ledger={attenuated},slack={slack}"
    )


if __name__ == "__main__":
    main()
