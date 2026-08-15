#!/usr/bin/env python3
"""Regressions for the new weighted double-endpoint fixed-split bounds."""

from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from alternating_double_endpoint_spectrum import (
    pair_arrays,
    scaled_endpoint_weights,
)
from mixed_endpoint_weighted_bound import (
    bound,
    improved_deterministic_ledger,
    mixed_fourier_value,
    refined_q32_certificate,
    sylvester,
    vertical_vertical_coefficient,
)
from weighted_same_orientation_certificate import certificate


def direct_mixed_q2(
    row_difference_law: np.ndarray,
    column_difference_law: np.ndarray,
) -> float:
    order = 2
    dimension = 4
    pair_left, pair_right = pair_arrays(dimension)
    pair_difference = pair_left ^ pair_right
    pair_count = len(pair_left)
    hadamard = sylvester(dimension) / 2
    endpoint = np.empty((dimension, dimension, pair_count))
    for singleton in range(dimension):
        endpoint[singleton] = (
            hadamard[singleton ^ pair_difference].T
            * scaled_endpoint_weights(
                order, singleton, pair_left, pair_right
            )
        )

    # Rows (i,b,F), columns (E,c,d).
    kernel = np.einsum(
        "ibe,bc,dcf->ibfecd",
        endpoint,
        hadamard,
        endpoint,
    ).reshape(
        dimension * dimension * pair_count,
        pair_count * dimension * dimension,
    )
    row_law = np.array(
        [
            2
            * row_difference_law[int(pair_difference[pair]) - 1]
            / dimension**3
            for singleton in range(dimension)
            for middle in range(dimension)
            for pair in range(pair_count)
        ]
    )
    column_law = np.array(
        [
            2
            * column_difference_law[int(pair_difference[pair]) - 1]
            / dimension**3
            for pair in range(pair_count)
            for middle in range(dimension)
            for singleton in range(dimension)
        ]
    )
    weighted = (
        np.sqrt(row_law)[:, None]
        * kernel
        * np.sqrt(column_law)[None, :]
    )
    return float(np.linalg.svd(weighted, compute_uv=False).sum())


def check_q2_reduction() -> None:
    random = np.random.default_rng(20260714)
    for _ in range(4):
        row_law = random.dirichlet(np.ones(3))
        column_law = random.dirichlet(np.ones(3))
        direct = direct_mixed_q2(row_law, column_law)
        reduced = mixed_fourier_value(2, row_law, column_law)
        if not np.isclose(direct, reduced, atol=3e-14):
            raise AssertionError(("mixed q2 reduction", direct, reduced))


def check_vertical_formula() -> None:
    for order in (4, 8):
        dimension = order * order
        vertical = np.arange(order, dimension, order) - 1
        law = np.zeros(dimension - 1)
        law[vertical] = 1 / (order - 1)
        direct = mixed_fourier_value(order, law, law)
        exact = float(vertical_vertical_coefficient(order))
        if not np.isclose(direct, exact, atol=4e-13):
            raise AssertionError(
                ("vertical mixed formula", order, direct, exact)
            )


def main() -> None:
    check_q2_reduction()
    check_vertical_formula()

    same = certificate(32)
    if same.supporting_upper >= 0.010905:
        raise AssertionError(
            ("same-orientation q32 upper", same.supporting_upper)
        )

    mixed = bound(32)
    if mixed.row_sum_upper >= mixed.old_fixed_split_bound * Decimal("0.25"):
        raise AssertionError(
            (
                "mixed q32 improvement",
                mixed.row_sum_upper,
                mixed.old_fixed_split_bound,
            )
        )
    required_local_threshold = (
        mixed.old_fixed_split_bound
        * Decimal("0.280707949532")
    )
    if mixed.row_sum_upper >= required_local_threshold:
        raise AssertionError(
            ("mixed q32 local threshold", mixed.row_sum_upper)
        )
    refined = refined_q32_certificate()
    if refined.advertised_upper >= Fraction(20344, 10**6):
        raise AssertionError(("refined mixed q32", refined))

    raw, attenuated, margin = improved_deterministic_ledger()
    if not 0.096 < margin - attenuated < 0.098:
        raise AssertionError(
            ("updated fixed-split slack", raw, attenuated, margin)
        )

    print(
        "weighted double-endpoint contraction passed: "
        f"same_q32={same.supporting_upper:.12g},"
        f"mixed_q32={float(refined.advertised_upper):.12g},"
        f"updated_attenuated_ledger={attenuated:.12g},"
        f"remaining_slack={margin-attenuated:.12g}"
    )


if __name__ == "__main__":
    main()
