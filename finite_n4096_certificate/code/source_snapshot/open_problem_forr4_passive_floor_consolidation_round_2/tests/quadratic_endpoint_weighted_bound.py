#!/usr/bin/env python3
"""Finite checks of the weighted quadratic-bent M_{a,1} contraction."""

from __future__ import annotations

from itertools import combinations

import numpy as np


SEED = 2026071420
N = 16


def sylvester(order: int) -> np.ndarray:
    result = np.array([[1]], dtype=int)
    while result.shape[0] < order:
        result = np.block([[result, result], [result, -result]])
    return result


def all_bent_pairs() -> tuple[np.ndarray, np.ndarray]:
    hadamard = sylvester(N)
    coordinates = np.arange(N)
    x_rows = []
    y_rows = []
    for mask in range(1 << N):
        x = (1 - 2 * ((mask >> coordinates) & 1)).astype(np.int8)
        transform = hadamard @ x
        if np.all(np.abs(transform) == 4):
            x_rows.append(x)
            y_rows.append((transform // 4).astype(np.int8))
    return np.array(x_rows), np.array(y_rows)


def main() -> None:
    rng = np.random.default_rng(SEED)
    x_values, y_values = all_bent_pairs()
    worst_ratio = 0.0
    checked = 0

    for degree in (1, 3, 5, 7, 9, 11):
        supports = list(combinations(range(N), degree))
        x_features = np.empty((len(x_values), len(supports)), dtype=np.int8)
        for column, support in enumerate(supports):
            x_features[:, column] = np.prod(x_values[:, support], axis=1)
        moment = x_features.T.astype(float) @ y_values / len(x_values)

        # Verify the exact XOR-row form before testing weighted norms.
        hadamard = sylvester(N) / 4
        for row, support in enumerate(supports):
            xor = 0
            for coordinate in support:
                xor ^= coordinate
            nonzero = np.flatnonzero(np.abs(hadamard[xor]) > 0)
            ratio = moment[row, nonzero] / hadamard[xor, nonzero]
            if not np.allclose(ratio, ratio[0], atol=2e-12):
                raise AssertionError(("XOR-row form", degree, support))
            if abs(ratio[0]) > 1 + 2e-12:
                raise AssertionError(("row multiplier", degree, support, ratio[0]))

        for _ in range(12):
            row_mass = float(rng.random())
            column_mass = float(rng.random())
            p = row_mass * rng.dirichlet(np.ones(len(supports)))
            q = column_mass * rng.dirichlet(np.ones(N))
            weighted = np.sqrt(p)[:, None] * moment * np.sqrt(q)[None, :]
            nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
            bound = np.sqrt(row_mass * column_mass)
            if nuclear > bound * (1 + 3e-11):
                raise AssertionError(
                    ("weighted endpoint bound", degree, nuclear, bound)
                )
            if bound:
                worst_ratio = max(worst_ratio, nuclear / bound)
            checked += 1

    print(
        "quadratic endpoint weighted bound passed: "
        f"bent_functions={len(x_values)}, weighted_instances={checked}, "
        f"worst_ratio={worst_ratio:.12g}"
    )


if __name__ == "__main__":
    main()
