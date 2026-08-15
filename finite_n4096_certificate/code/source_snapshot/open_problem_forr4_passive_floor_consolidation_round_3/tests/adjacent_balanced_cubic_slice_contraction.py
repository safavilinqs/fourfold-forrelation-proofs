#!/usr/bin/env python3
"""Regression for the second balanced chain-aware contraction."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from adjacent_balanced_cubic_slice_contraction import (  # noqa: E402
    adjacent_balanced_coefficient,
    adjacent_balanced_contraction,
    adjacent_balanced_orbit_entries,
)
from adjacent_cubic_quintic_orbit_witness import (  # noqa: E402
    exact_link_moments,
)
from leading_balanced_disjointness_contraction import (  # noqa: E402
    disjointness_schur_factor,
)
from occupation_compatible_sector_optimization import (  # noqa: E402
    endpoint_singleton_slice_energies,
)


SEED = 2026071502


def cyclic_fourier(dimension: int) -> np.ndarray:
    indices = np.arange(dimension)
    return np.exp(
        2j * np.pi * indices[:, None] * indices[None, :] / dimension
    ) / np.sqrt(dimension)


def direct_abstract_chain_checks() -> float:
    """Stress the Schur composition on complete small tensors."""

    rng = np.random.default_rng(SEED)
    dimension = 6
    pairs = tuple(combinations(range(dimension), 2))
    fours = tuple(combinations(range(dimension), 4))
    quintic_mask = np.asarray(
        [
            [float(element not in support) for element in range(dimension)]
            for support in fours
        ]
    )
    fourier = cyclic_fourier(dimension)
    slice_energy = 0.23
    bound = (
        disjointness_schur_factor(dimension, 4)
        * np.sqrt(slice_energy)
        / np.sqrt(dimension)
    )
    worst = 0.0

    for _ in range(5):
        # One column of the occurrence symbol is a vector over x.  Normalize
        # every (b,E) vector to the theorem's squared-slice cap.
        cubic = np.zeros((dimension, dimension, len(pairs)), dtype=complex)
        for b in range(dimension):
            for pair_index, pair in enumerate(pairs):
                values = (
                    rng.normal(size=dimension)
                    + 1j * rng.normal(size=dimension)
                )
                values[list(pair)] = 0
                values *= np.sqrt(slice_energy) / np.linalg.norm(values)
                cubic[:, b, pair_index] = values

        latent = 7
        law = rng.dirichlet(np.ones(latent))
        row_features = np.exp(
            2j
            * np.pi
            * rng.random((latent, dimension, len(fours)))
        )
        column_features = np.exp(
            2j
            * np.pi
            * rng.random((latent, len(pairs), dimension))
        )
        completed_adjacent = np.einsum(
            "lxf,l,lpe->xfpe",
            row_features.conj(),
            law,
            column_features,
        )
        adjacent = completed_adjacent * quintic_mask[None, :, None, :]

        # Rows are (x,F), columns are (a,b,E,e).  The central Fourier entry
        # is a column-only scalar of modulus 1/sqrt(N).
        tensor = np.einsum(
            "ab,xbp,xfpe->xfabpe",
            fourier,
            cubic,
            adjacent,
        ).reshape(
            dimension * len(fours),
            dimension * dimension * len(pairs) * dimension,
        )
        row_law = rng.dirichlet(np.ones(tensor.shape[0]))
        column_law = rng.dirichlet(np.ones(tensor.shape[1]))
        weighted = (
            np.sqrt(row_law)[:, None]
            * tensor
            * np.sqrt(column_law)[None, :]
        )
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > bound * (1 + 4e-12):
            raise AssertionError(("abstract adjacent contraction", nuclear, bound))
        worst = max(worst, nuclear / bound)
    return worst


def main() -> None:
    # Complete q=4 signed-permutation enumeration independently checks that
    # the fixed-pair M_13 slice is the exact cubic table entry T_2.
    moments = exact_link_moments(4)
    support_index = {
        support: index
        for index, support in enumerate(moments.supports_three)
    }
    dimension = 16
    maximum = 0.0
    for singleton in range(dimension):
        for pair in combinations(range(dimension), 2):
            extensions = [
                support_index[tuple(sorted(pair + (cell,)))]
                for cell in range(dimension)
                if cell not in pair
            ]
            maximum = max(
                maximum,
                float(np.sum(moments.moment_13[singleton, extensions] ** 2)),
            )
    expected = endpoint_singleton_slice_energies(4)[2]
    if not np.isclose(maximum, expected, atol=2e-14):
        raise AssertionError(("exact cubic fixed-pair slice", maximum, expected))

    orbit = adjacent_balanced_orbit_entries()
    if len(orbit) != 4:
        raise AssertionError(("adjacent balanced orbit", orbit))
    coefficient = adjacent_balanced_coefficient()
    if not np.isclose(coefficient, 0.016272469279635413, atol=2e-15):
        raise AssertionError(("adjacent balanced coefficient", coefficient))
    result = adjacent_balanced_contraction()
    if not np.isclose(result.optimized_total, 0.32556385796981924, atol=3e-11):
        raise AssertionError(("updated adjacent ledger", result))
    if result.threshold_slack <= 0.00776:
        raise AssertionError(("updated adjacent slack", result))
    if not result.coefficient < result.provisional_coefficient:
        raise AssertionError(("coefficient does not beat 1/q", result))
    if not result.generic_two_mask_coefficient > 0.225:
        raise AssertionError(("generic obstruction missing", result))

    worst = direct_abstract_chain_checks()
    print(
        "adjacent balanced cubic-slice contraction passed: "
        f"cubic_slice={result.cubic_pair_slice_energy:.12g},"
        f"coefficient={result.coefficient:.12g},"
        f"generic_two_mask={result.generic_two_mask_coefficient:.12g},"
        f"ledger_total={result.optimized_total:.12g},"
        f"threshold_slack={result.threshold_slack:.12g},"
        f"worst_random_ratio={worst:.12g}"
    )


if __name__ == "__main__":
    main()
