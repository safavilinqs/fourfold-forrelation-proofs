#!/usr/bin/env python3
"""Regression for the leading balanced distinct-label contraction."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from leading_balanced_disjointness_contraction import (  # noqa: E402
    disjointness_schur_factor,
    leading_balanced_coefficient,
    leading_balanced_contraction,
    leading_balanced_orbit_entries,
)


SEED = 2026071501


def disjointness_factors(
    dimension: int,
    selected: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Construct the centered factorization used by the theorem."""

    subsets = tuple(combinations(range(dimension), selected))
    alpha = 1 - selected / dimension
    row_energy = selected * alpha
    column_energy = 1 - 1 / dimension
    scale = (column_energy / row_energy) ** 0.25
    rows = np.empty((len(subsets), dimension + 1))
    columns = np.empty((dimension, dimension + 1))
    rows[:, 0] = np.sqrt(alpha)
    columns[:, 0] = np.sqrt(alpha)
    for index, subset in enumerate(subsets):
        centered = np.full(dimension, -selected / dimension)
        centered[list(subset)] += 1
        rows[index, 1:] = scale * centered
    columns[:, 1:] = -(np.eye(dimension) - 1 / dimension) / scale
    return rows, columns


def cyclic_fourier(dimension: int) -> np.ndarray:
    indices = np.arange(dimension)
    return np.exp(
        2j * np.pi * indices[:, None] * indices[None, :] / dimension
    ) / np.sqrt(dimension)


def direct_random_chain_checks() -> float:
    """Stress the abstract theorem on a fully constructed small tensor."""

    rng = np.random.default_rng(SEED)
    # Dimension six keeps the full matrix small while making the 4|1
    # disjointness Schur factor strictly larger than one.
    dimension = 6
    selected = 4
    triples = tuple(combinations(range(dimension), 3))
    fours = tuple(combinations(range(dimension), selected))
    four_masks = np.asarray(
        [
            [float(element not in support) for element in range(dimension)]
            for support in fours
        ]
    )
    fourier = cyclic_fourier(dimension)
    schur = disjointness_schur_factor(dimension, selected)
    worst = 0.0

    for _ in range(8):
        cubic_weights = rng.uniform(-1, 1, len(triples))
        cubic_labels = np.asarray(
            [sum(support) % dimension for support in triples]
        )
        cubic = (
            cubic_weights[:, None]
            * fourier[cubic_labels]
        )

        latent = 7
        latent_law = rng.dirichlet(np.ones(latent))
        row_features = np.exp(
            2j * np.pi * rng.random((latent, len(fours)))
        )
        column_features = np.exp(
            2j * np.pi
            * rng.random((latent, dimension, dimension))
        )
        completed_quintic = np.einsum(
            "lf,l,lec->fec",
            row_features.conj(),
            latent_law,
            column_features,
        )
        quintic = completed_quintic * four_masks[:, :, None]

        # Rows are (b,F), columns are (Q,c,e).
        tensor = np.einsum(
            "qb,bc,fec->bfqce",
            cubic,
            fourier,
            quintic,
        ).reshape(
            dimension * len(fours),
            len(triples) * dimension * dimension,
        )
        row_law = rng.dirichlet(np.ones(tensor.shape[0]))
        column_law = rng.dirichlet(np.ones(tensor.shape[1]))
        weighted = (
            np.sqrt(row_law)[:, None]
            * tensor
            * np.sqrt(column_law)[None, :]
        )
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        bound = schur / np.sqrt(dimension)
        if nuclear > bound * (1 + 3e-12):
            raise AssertionError(("random chain contraction", nuclear, bound))
        worst = max(worst, nuclear / bound)
    return worst


def main() -> None:
    for dimension, selected in ((5, 1), (5, 4), (8, 2), (8, 4)):
        rows, columns = disjointness_factors(dimension, selected)
        subsets = tuple(combinations(range(dimension), selected))
        expected = np.asarray(
            [
                [float(element not in subset) for element in range(dimension)]
                for subset in subsets
            ]
        )
        if not np.allclose(rows @ columns.T, expected, atol=2e-14):
            raise AssertionError(("disjointness factorization", dimension, selected))
        factor = max(np.linalg.norm(rows, axis=1)) * max(
            np.linalg.norm(columns, axis=1)
        )
        expected_factor = disjointness_schur_factor(dimension, selected)
        if not np.isclose(factor, expected_factor, atol=2e-14):
            raise AssertionError(
                ("disjointness factor norm", dimension, selected, factor)
            )

    orbit = leading_balanced_orbit_entries()
    if len(orbit) != 4:
        raise AssertionError(("leading orbit", orbit))
    coefficient = leading_balanced_coefficient()
    if not np.isclose(coefficient, 0.09347527457750368, atol=2e-15):
        raise AssertionError(("leading coefficient", coefficient))
    result = leading_balanced_contraction()
    if not np.isclose(result.optimized_total, 0.32628186091832256, atol=3e-11):
        raise AssertionError(("updated ledger total", result))
    if result.threshold_slack <= 0.00705:
        raise AssertionError(("updated ledger slack", result))

    worst = direct_random_chain_checks()
    print(
        "leading balanced disjointness contraction passed: "
        f"coefficient={result.coefficient:.12g},"
        f"ledger_total={result.optimized_total:.12g},"
        f"threshold_slack={result.threshold_slack:.12g},"
        f"worst_random_ratio={worst:.12g}"
    )


if __name__ == "__main__":
    main()
