#!/usr/bin/env python3
"""Stress the complete weighted cut table for a three-link path tensor."""

from __future__ import annotations

from itertools import product

import numpy as np


SEED = 2026071422


def moment_matrix(
    rng: np.random.Generator, latent: int, left: int, right: int
) -> np.ndarray:
    probability = rng.dirichlet(np.ones(latent))
    left_features = np.exp(2j * np.pi * rng.random((latent, left)))
    right_features = np.exp(2j * np.pi * rng.random((latent, right)))
    return left_features.conj().T @ (probability[:, None] * right_features)


def flatten(tensor: np.ndarray, mask: int) -> np.ndarray:
    row_axes = [axis for axis in range(4) if (mask >> axis) & 1]
    column_axes = [axis for axis in range(4) if not ((mask >> axis) & 1)]
    ordered = np.transpose(tensor, row_axes + column_axes)
    row_size = int(np.prod([tensor.shape[axis] for axis in row_axes]))
    column_size = int(
        np.prod([tensor.shape[axis] for axis in column_axes])
    )
    return ordered.reshape(row_size, column_size)


def cut_bound(
    mask: int,
    kappas: tuple[float, float, float],
    operators: tuple[float, float, float],
    dimensions: tuple[int, int, int, int],
) -> float:
    if mask > (15 ^ mask):
        mask ^= 15
    k1, k2, k3 = kappas
    c1, _, c3 = operators
    _, d2, d3, _ = dimensions
    tau12 = min(1.0, np.sqrt(d2) * k1 * k2)
    tau23 = min(1.0, np.sqrt(d3) * k2 * k3)
    table = {
        0: k1 * k2 * k3,
        1: k2 * k3,
        2: k3 * tau12,
        3: k1 * k3,
        4: k1 * tau23,
        5: c1 * k2 * c3,
        6: k2,
        7: k1 * k2,
    }
    return table[mask]


def main() -> None:
    rng = np.random.default_rng(SEED)
    worst_ratio = 0.0
    worst_alternating_coherence_ratio = 0.0
    checked = 0
    for dimensions in ((2, 3, 2, 4), (3, 2, 4, 2), (4, 3, 3, 2)):
        d1, d2, d3, d4 = dimensions
        for _ in range(30):
            matrices = (
                moment_matrix(rng, 7, d1, d2),
                moment_matrix(rng, 8, d2, d3),
                moment_matrix(rng, 9, d3, d4),
            )
            m1, m2, m3 = matrices
            tensor = (
                m1[:, :, None, None]
                * m2[None, :, :, None]
                * m3[None, None, :, :]
            )
            kappas = tuple(float(np.max(np.abs(matrix))) for matrix in matrices)
            operators = tuple(
                float(np.linalg.svd(matrix, compute_uv=False)[0])
                for matrix in matrices
            )
            outer_row_energy = float(np.max(np.sum(np.abs(m1) ** 2, axis=1)))
            outer_column_energy = float(
                np.max(np.sum(np.abs(m3) ** 2, axis=0))
            )
            for mask in range(16):
                matrix = flatten(tensor, mask)
                row_mass = float(rng.random())
                column_mass = float(rng.random())
                p = row_mass * rng.dirichlet(np.ones(matrix.shape[0]))
                q = column_mass * rng.dirichlet(np.ones(matrix.shape[1]))
                weighted = np.sqrt(p)[:, None] * matrix * np.sqrt(q)[None, :]
                nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
                bound = cut_bound(mask, kappas, operators, dimensions)
                bound *= np.sqrt(row_mass * column_mass)
                if nuclear > bound * (1 + 2e-11):
                    raise AssertionError(
                        (
                            "three-link cut bound",
                            dimensions,
                            mask,
                            nuclear,
                            bound,
                        )
                    )
                if bound:
                    worst_ratio = max(worst_ratio, nuclear / bound)
                if mask in (5, 10):
                    coherence_bound = kappas[1] * np.sqrt(
                        row_mass * column_mass
                    )
                    if coherence_bound:
                        worst_alternating_coherence_ratio = max(
                            worst_alternating_coherence_ratio,
                            nuclear / coherence_bound,
                        )
                    bessel_bound = (
                        np.sqrt(outer_row_energy * outer_column_energy)
                        * kappas[1]
                        * np.sqrt(row_mass * column_mass)
                    )
                    if nuclear > bessel_bound * (1 + 2e-11):
                        raise AssertionError(
                            (
                                "Bessel-refined alternating cut",
                                dimensions,
                                mask,
                                nuclear,
                                bessel_bound,
                            )
                        )
                checked += 1

    # The singleton Hadamard sector must reproduce the exact graph-cut
    # powers used in the one-batch minimal-chain certificate.
    for n in (2, 4, 8, 16):
        kappa = n ** -0.5
        kappas = (kappa, kappa, kappa)
        operators = (1.0, 1.0, 1.0)
        dimensions = (n, n, n, n)
        expected = {}
        for mask in range(16):
            representative = min(mask, 15 ^ mask)
            if representative == 0:
                expected[mask] = n ** -1.5
            elif representative in (5, 6):
                expected[mask] = n ** -0.5
            else:
                expected[mask] = n ** -1.0
            observed = cut_bound(mask, kappas, operators, dimensions)
            if not np.isclose(observed, expected[mask], atol=2e-15):
                raise AssertionError(
                    ("Hadamard cut table", n, mask, observed, expected[mask])
                )

    print(
        "three-link weighted path contraction passed: "
        f"random_cuts={checked}, worst_ratio={worst_ratio:.12g}, "
        "worst_alternating_coherence_ratio="
        f"{worst_alternating_coherence_ratio:.12g}, "
        "Hadamard_dimensions=2,4,8,16"
    )


if __name__ == "__main__":
    main()
