#!/usr/bin/env python3
"""Vary the invariant endpoint-difference correlation in the physical law.

The Round 3 lower witness sampled the two nonzero vertical pair differences
independently and uniformly.  Under the full row-label linear symmetry, an
ordered pair of nonzero differences has two orbits: equal and distinct.
This search varies the total mass on the equal orbit and evaluates the same
exact fourteen-class Fourier reduction.

The class reduction is valid only for these symmetry-invariant laws.  This
module deliberately exposes no arbitrary pair-law evaluator: nonsymmetric
laws require the unreduced frequency sum.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import sqrt
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parents[1]
ROUND3 = ROOT.parent / "open_problem_forr4_passive_floor_consolidation_round_3"
sys.path.insert(0, str(ROUND3 / "searches"))

from opposite_endpoint_orbit_scan import (  # noqa: E402
    cubic_response,
    quintic_response,
    support_xor,
    triple_orbit_representatives,
    walsh_transform,
)
from opposite_endpoint_vertical_mixture_witness import character  # noqa: E402


@dataclass(frozen=True)
class CorrelatedVerticalMixture:
    order: int
    equal_difference_mass: float
    distinct_difference_mass: float
    triple_orbits: int
    frequency_classes: int
    coefficient: float


@dataclass(frozen=True)
class CorrelatedMixtureSearch:
    order: int
    independent_equal_mass: float
    independent_coefficient: float
    best_equal_mass: float
    best_coefficient: float
    improvement: float


@lru_cache(maxsize=None)
def _spectral_blocks(
    order: int,
) -> tuple[tuple[tuple[np.ndarray, int], ...], int]:
    """Build the fourteen unweighted invariant Fourier blocks."""

    if order < 4 or order & (order - 1):
        raise ValueError(("power-of-two order at least four required", order))
    q = order
    dimension = q * q
    vertical = tuple(row * q for row in range(1, q))
    triples = tuple(
        tuple(row * q for row in triple)
        for triple in triple_orbit_representatives(q)
    )
    cubic = np.asarray(
        [
            walsh_transform(cubic_response(q, difference, False))
            for difference in vertical
        ]
    )
    quintic = np.empty((q - 1, len(triples), dimension))
    frequencies = range(dimension)
    for difference_index, difference in enumerate(vertical):
        for triple_index, triple in enumerate(triples):
            triple_xor = support_xor(triple)
            quintic[difference_index, triple_index] = walsh_transform(
                quintic_response(q, difference, triple, False)
            ) * np.asarray(
                [character(triple_xor, frequency) for frequency in frequencies]
            )

    row_frequency_classes = (
        (0, 0, 1),
        (0, 1, q - 1),
        (1, 0, q - 1),
        (1, 2, (q - 1) * (q // 2 - 1)),
        (1, 1, (q - 1) * (q // 2)),
    )
    blocks = []
    for alpha_row, gamma_row, row_multiplicity in row_frequency_classes:
        alpha_columns = (
            ((0, 1), (1, q - 1))
            if alpha_row == 0
            else ((0, q),)
        )
        for alpha_column, alpha_multiplicity in alpha_columns:
            for gamma_column, gamma_multiplicity in ((0, 1), (1, q - 1)):
                alpha = alpha_row * q + alpha_column
                gamma = gamma_row * q + gamma_column
                shifted_frequencies = np.asarray(
                    [gamma ^ difference for difference in vertical]
                )
                endpoint = quintic[:, :, shifted_frequencies].transpose(2, 0, 1)
                matrix = (
                    cubic[:, alpha, None, None] * endpoint
                ).reshape((q - 1) ** 2, len(triples))
                multiplicity = (
                    row_multiplicity
                    * alpha_multiplicity
                    * gamma_multiplicity
                )
                blocks.append((matrix, multiplicity))
    return tuple(blocks), len(triples)


def correlated_vertical_mixture(
    order: int,
    equal_difference_mass: float,
) -> CorrelatedVerticalMixture:
    """Evaluate one full-row-symmetry-invariant correlated physical law."""

    if not 0 <= equal_difference_mass <= 1:
        raise ValueError(("equal mass outside probability simplex", equal_difference_mass))
    q = order
    dimension = q * q
    blocks, triple_count = _spectral_blocks(order)
    pair_probabilities = np.full(
        (q - 1, q - 1),
        (1 - equal_difference_mass) / ((q - 1) * (q - 2)),
    )
    np.fill_diagonal(
        pair_probabilities,
        equal_difference_mass / (q - 1),
    )
    if not np.isclose(pair_probabilities.sum(), 1, atol=2e-15):
        raise AssertionError(("row law normalization", pair_probabilities.sum()))
    row_amplitudes = np.sqrt(pair_probabilities).reshape(-1)
    column_amplitude = 1 / sqrt(triple_count)
    nuclear_sum = 0.0
    for matrix, multiplicity in blocks:
        weighted = matrix * row_amplitudes[:, None] * column_amplitude
        nuclear_sum += multiplicity * float(
            np.linalg.svd(weighted, compute_uv=False).sum()
        )
    return CorrelatedVerticalMixture(
        order=q,
        equal_difference_mass=equal_difference_mass,
        distinct_difference_mass=1 - equal_difference_mass,
        triple_orbits=triple_count,
        frequency_classes=len(blocks),
        coefficient=nuclear_sum / dimension**3,
    )


def search_correlated_vertical_mixture(
    order: int = 32,
) -> CorrelatedMixtureSearch:
    """Numerically maximize the one-parameter invariant physical family."""

    independent_mass = 1 / (order - 1)
    independent = correlated_vertical_mixture(order, independent_mass)
    cache: dict[float, float] = {independent_mass: independent.coefficient}

    def coefficient(mass: float) -> float:
        key = float(mass)
        if key not in cache:
            cache[key] = correlated_vertical_mixture(order, key).coefficient
        return cache[key]

    grid = np.linspace(0, 1, 41)
    best_index = int(np.argmax([coefficient(float(mass)) for mass in grid]))
    left = float(grid[max(0, best_index - 1)])
    right = float(grid[min(len(grid) - 1, best_index + 1)])
    optimum = minimize_scalar(
        lambda mass: -coefficient(float(mass)),
        bounds=(left, right),
        method="bounded",
        options={"xatol": 1e-13},
    )
    candidates = set(cache) | {0.0, 1.0, float(optimum.x)}
    best_mass = max(candidates, key=coefficient)
    best = coefficient(best_mass)
    return CorrelatedMixtureSearch(
        order=order,
        independent_equal_mass=independent_mass,
        independent_coefficient=independent.coefficient,
        best_equal_mass=best_mass,
        best_coefficient=best,
        improvement=best - independent.coefficient,
    )


def main() -> None:
    result = search_correlated_vertical_mixture()
    print(
        "opposite-endpoint invariant correlated mixture: "
        f"q={result.order},"
        f"independent_mass={result.independent_equal_mass:.15g},"
        f"independent={result.independent_coefficient:.15g},"
        f"best_mass={result.best_equal_mass:.15g},"
        f"best={result.best_coefficient:.15g},"
        f"improvement={result.improvement:.15g}"
    )


if __name__ == "__main__":
    main()
