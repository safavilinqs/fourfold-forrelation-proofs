#!/usr/bin/env python3
"""Random and saturating checks of the Gram-dressed two-link tail bound."""

from __future__ import annotations

import numpy as np


SEED = 2026071423


def sylvester(order: int) -> np.ndarray:
    result = np.array([[1.0]])
    while result.shape[0] < order:
        result = np.block([[result, result], [result, -result]])
    return result / np.sqrt(order)


def main() -> None:
    rng = np.random.default_rng(SEED)
    worst_ratio = 0.0
    checked = 0
    for dimension in (2, 4, 8, 16):
        hadamard = sylvester(dimension)
        wedge = np.einsum("bc,cd->bdc", hadamard, hadamard).reshape(
            dimension * dimension, dimension
        )
        target = dimension ** -0.5
        for row_repetitions, column_repetitions in ((1, 1), (2, 3), (4, 2)):
            base_rows = np.repeat(
                np.arange(dimension * dimension), row_repetitions
            )
            base_columns = np.repeat(
                np.arange(dimension), column_repetitions
            )
            lifted = wedge[np.ix_(base_rows, base_columns)]
            for _ in range(12):
                latent = 9
                probability = rng.dirichlet(np.ones(latent))
                left_features = np.exp(
                    2j * np.pi * rng.random((latent, lifted.shape[0]))
                )
                right_features = np.exp(
                    2j * np.pi * rng.random((latent, lifted.shape[1]))
                )
                gram = (
                    left_features.conj().T
                    @ (probability[:, None] * right_features)
                )
                dressed = gram * lifted
                row_mass = float(rng.random())
                column_mass = float(rng.random())
                p = row_mass * rng.dirichlet(np.ones(lifted.shape[0]))
                q = column_mass * rng.dirichlet(np.ones(lifted.shape[1]))
                weighted = (
                    np.sqrt(p)[:, None] * dressed * np.sqrt(q)[None, :]
                )
                nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
                bound = target * np.sqrt(row_mass * column_mass)
                if nuclear > bound * (1 + 3e-11):
                    raise AssertionError(
                        (
                            "Gram-dressed tail contraction",
                            dimension,
                            row_repetitions,
                            column_repetitions,
                            nuclear,
                            bound,
                        )
                    )
                if bound:
                    worst_ratio = max(worst_ratio, nuclear / bound)
                checked += 1

        # With no dressing and uniform base weights, the wedge has N
        # singular values one and saturates the coefficient 1/sqrt(N).
        p = np.full(dimension * dimension, 1 / (dimension * dimension))
        q = np.full(dimension, 1 / dimension)
        weighted = np.sqrt(p)[:, None] * wedge * np.sqrt(q)[None, :]
        saturation = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if not np.isclose(saturation, target, atol=3e-14):
            raise AssertionError(
                ("two-link wedge saturation", dimension, saturation, target)
            )

    print(
        "Gram-dressed tail contraction passed: "
        f"instances={checked}, worst_random_ratio={worst_ratio:.12g}, "
        "saturating_dimensions=2,4,8,16"
    )


if __name__ == "__main__":
    main()
