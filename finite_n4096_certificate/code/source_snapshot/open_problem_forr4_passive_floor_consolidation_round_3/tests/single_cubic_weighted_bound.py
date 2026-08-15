#!/usr/bin/env python3
"""Regression for the endpoint single-cubic weighted coefficients."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from alternating_double_endpoint_spectrum import scaled_endpoint_weights
from mixed_endpoint_weighted_bound import sylvester
from single_cubic_weighted_bound import bound


def direct_q2_values() -> dict[int, float]:
    order = 2
    dimension = 4
    hadamard = sylvester(dimension) / 2
    pairs = list(combinations(range(dimension), 2))
    left = np.array([pair[0] for pair in pairs], dtype=np.int32)
    right = np.array([pair[1] for pair in pairs], dtype=np.int32)
    difference = left ^ right
    endpoint = np.empty((dimension, len(pairs), dimension))
    for singleton in range(dimension):
        endpoint[singleton] = (
            scaled_endpoint_weights(order, singleton, left, right)[:, None]
            * hadamard[singleton ^ difference]
        )
    tensor = np.einsum(
        "ieb,bc,cd->iebcd", endpoint, hadamard, hadamard
    )
    values: dict[int, float] = {}
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
            int(np.prod([tensor.shape[axis] for axis in row_axes])),
            int(np.prod([tensor.shape[axis] for axis in column_axes])),
        )
        values[mask] = float(
            np.linalg.svd(matrix, compute_uv=False).sum()
            / np.sqrt(matrix.shape[0] * matrix.shape[1])
        )
    return values


def main() -> None:
    q2 = bound(2)
    values = direct_q2_values()
    extreme = max(
        value for mask, value in values.items() if mask.bit_count() in (0, 3)
    )
    balanced = max(
        value for mask, value in values.items() if mask.bit_count() in (1, 2)
    )
    if not np.isclose(
        extreme, float(q2.extreme_singletons), atol=3e-14
    ):
        raise AssertionError(("q2 extreme cubic", extreme, q2))
    if not np.isclose(
        balanced, float(q2.balanced_singletons), atol=3e-14
    ):
        raise AssertionError(("q2 balanced cubic", balanced, q2))

    q32 = bound(32)
    if not 0.000324 < float(q32.extreme_singletons) < 0.000326:
        raise AssertionError(("q32 extreme cubic", q32))
    if not 0.01039 < float(q32.balanced_singletons) < 0.01040:
        raise AssertionError(("q32 balanced cubic", q32))
    print(
        "single-cubic weighted bound passed: "
        f"extreme_q32={q32.extreme_singletons},"
        f"balanced_q32={q32.balanced_singletons}"
    )


if __name__ == "__main__":
    main()
