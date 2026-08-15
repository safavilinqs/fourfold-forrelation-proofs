#!/usr/bin/env python3
"""Regression for the fourth balanced shared-law contraction."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from internal_singleton_shared_law_contraction import (  # noqa: E402
    cubic_twirled_coefficient,
    internal_singleton_coefficient,
    internal_singleton_contraction,
    internal_singleton_orbit_entries,
    vertical_mixture_diagnostic,
)
from opposite_endpoint_orbit_scan import (  # noqa: E402
    cubic_weight,
    endpoint_moment,
    quintic_weight,
    support_xor,
    sylvester,
)


SEED = 2026071504


def pair_type(pair: tuple[int, int], order: int) -> int:
    first_row, first_column = divmod(pair[0], order)
    second_row, second_column = divmod(pair[1], order)
    if first_column == second_column:
        return 0
    if first_row == second_row:
        return 1
    return 2


def exact_q4_cubic_twirl() -> tuple[float, float]:
    """Build the complete q=4 cubic matrix and its twirled optimizer."""

    q = 4
    dimension = q * q
    pairs = tuple(combinations(range(dimension), 2))
    hadamard = sylvester(dimension)
    matrix = np.zeros((dimension * dimension, len(pairs)))
    for x in range(dimension):
        for b in range(dimension):
            row = x * dimension + b
            for pair_index, pair in enumerate(pairs):
                if x in pair:
                    continue
                support = tuple(sorted(pair + (x,)))
                matrix[row, pair_index] = (
                    cubic_weight(support, q)
                    * hadamard[support_xor(support), b]
                )

    result = cubic_twirled_coefficient(q)
    contributions = np.asarray(
        [
            result.vertical_contribution,
            result.horizontal_contribution,
            result.general_contribution,
        ]
    )
    type_weights = contributions**2 / np.sum(contributions**2)
    type_counts = np.bincount(
        [pair_type(pair, q) for pair in pairs], minlength=3
    )
    row_law = np.full(matrix.shape[0], 1 / matrix.shape[0])
    column_law = np.asarray(
        [
            type_weights[pair_type(pair, q)]
            / type_counts[pair_type(pair, q)]
            for pair in pairs
        ]
    )
    weighted = (
        np.sqrt(row_law)[:, None]
        * matrix
        * np.sqrt(column_law)[None, :]
    )
    attained = float(np.linalg.svd(weighted, compute_uv=False).sum())

    rng = np.random.default_rng(SEED)
    worst_random = 0.0
    for _ in range(12):
        row_law = rng.dirichlet(np.ones(matrix.shape[0]))
        column_law = rng.dirichlet(np.ones(matrix.shape[1]))
        weighted = (
            np.sqrt(row_law)[:, None]
            * matrix
            * np.sqrt(column_law)[None, :]
        )
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        worst_random = max(worst_random, nuclear / result.coefficient)
        if nuclear > result.coefficient * (1 + 4e-12):
            raise AssertionError(("q4 cubic twirl", nuclear, result))
    return attained, worst_random


def exact_q4_quintic_completion() -> float:
    """Check B=completed-overlap and the overlap row-energy bound."""

    q = 4
    dimension = q * q
    pairs = tuple(combinations(range(dimension), 2))
    triples = tuple(combinations(range(dimension), 3))
    hadamard = sylvester(dimension)
    maximum_collision_energy = 0.0
    for pair in pairs:
        pair_set = set(pair)
        for singleton in range(dimension):
            collision_energy = 0.0
            for triple in triples:
                overlap = pair_set.intersection(triple)
                union = tuple(sorted(pair + triple))
                if not overlap:
                    actual = (
                        quintic_weight(union, q)
                        * hadamard[support_xor(union), singleton]
                    )
                    completed = actual
                else:
                    symmetric_difference = tuple(
                        sorted(pair_set.symmetric_difference(triple))
                    )
                    if len(symmetric_difference) == 3:
                        completed = endpoint_moment(
                            symmetric_difference,
                            singleton,
                            q,
                            degree=3,
                            high_only=False,
                        )
                    elif len(symmetric_difference) == 1:
                        completed = hadamard[
                            symmetric_difference[0], singleton
                        ]
                    else:
                        raise AssertionError(symmetric_difference)
                    actual = 0.0
                    collision_energy += completed * completed
                if not np.isclose(
                    actual,
                    completed - (completed if overlap else 0.0),
                    atol=2e-14,
                ):
                    raise AssertionError((pair, triple, singleton))
            maximum_collision_energy = max(
                maximum_collision_energy, collision_energy
            )
    if maximum_collision_energy > 2 + 3e-14:
        raise AssertionError(
            ("quintic overlap row energy", maximum_collision_energy)
        )
    return maximum_collision_energy


def main() -> None:
    attained, worst_random = exact_q4_cubic_twirl()
    q4_cubic = cubic_twirled_coefficient(4)
    if not np.isclose(attained, q4_cubic.coefficient, atol=3e-14):
        raise AssertionError(("twirled optimizer", attained, q4_cubic))
    collision_energy = exact_q4_quintic_completion()

    orbit = internal_singleton_orbit_entries()
    if len(orbit) != 4:
        raise AssertionError(("internal singleton orbit", orbit))
    cubic = cubic_twirled_coefficient()
    if (
        cubic.vertical_spectrum_numerator,
        cubic.horizontal_spectrum_numerator,
        cubic.general_spectrum_numerator,
    ) != (59644, 1984, 1984):
        raise AssertionError(("q32 cubic spectra", cubic))
    if not np.isclose(cubic.coefficient, 0.3326532036394109, atol=2e-15):
        raise AssertionError(("q32 cubic twirl", cubic))
    coefficient = internal_singleton_coefficient()
    if not np.isclose(coefficient, 0.02509674611853515, atol=2e-15):
        raise AssertionError(("shared-law coefficient", coefficient))

    result = internal_singleton_contraction()
    if not result.coefficient < result.provisional_coefficient:
        raise AssertionError(("provisional improvement", result))
    if not result.coefficient < result.acceptance_gate:
        raise AssertionError(("acceptance gate", result))
    if not result.simple_slice_threshold_overshoot > 0.00385:
        raise AssertionError(("simple slice obstruction", result))
    if not result.threshold_slack > 0.00096:
        raise AssertionError(("four-theorem ledger", result))
    if result.next_unresolved_entries[0] != (
        (3, 1, 1, 5),
        (0, 0, 1, 4),
    ):
        raise AssertionError(("reranked next orbit", result))
    if not np.isclose(
        result.next_admissible_coefficient,
        0.0542506297760259,
        atol=4e-14,
    ):
        raise AssertionError(("next admissible coefficient", result))
    if not np.isclose(
        result.vertical_mixture_diagnostic,
        0.006656191979278528,
        atol=2e-15,
    ):
        raise AssertionError(("vertical mixture diagnostic", result))
    if not np.isclose(
        vertical_mixture_diagnostic(4),
        0.0390625,
        atol=2e-15,
    ):
        raise AssertionError("q4 vertical mixture")

    print(
        "internal singleton shared-law contraction passed: "
        f"q4_cubic={q4_cubic.coefficient:.12g},"
        f"q4_collision_energy={collision_energy:.12g},"
        f"q4_random_ratio={worst_random:.12g},"
        f"q32_cubic={cubic.coefficient:.12g},"
        f"coefficient={result.coefficient:.12g},"
        f"simple_overshoot={result.simple_slice_threshold_overshoot:.12g},"
        f"ledger_total={result.optimized_total:.12g},"
        f"threshold_slack={result.threshold_slack:.12g},"
        f"next_admissible={result.next_admissible_coefficient:.12g},"
        f"vertical_mixture={result.vertical_mixture_diagnostic:.12g}"
    )


if __name__ == "__main__":
    main()
