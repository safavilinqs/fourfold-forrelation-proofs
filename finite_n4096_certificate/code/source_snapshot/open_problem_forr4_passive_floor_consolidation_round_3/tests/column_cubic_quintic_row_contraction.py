#!/usr/bin/env python3
"""Regression for the fifth balanced row-energy contraction."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from column_cubic_quintic_row_contraction import (  # noqa: E402
    column_cubic_quintic_coefficient,
    column_cubic_quintic_contraction,
    column_cubic_quintic_orbit_entries,
)
from occupation_compatible_sector_optimization import (  # noqa: E402
    endpoint_quintic_singleton_slice_energies,
)
from opposite_endpoint_orbit_scan import (  # noqa: E402
    cubic_weight,
    quintic_weight,
    support_xor,
    sylvester,
)


SEED = 2026071505


def exact_q4_quintic_row_energy() -> float:
    """Enumerate every fixed-four M_15 row directly from the moments."""

    order = 4
    dimension = order * order
    hadamard = sylvester(dimension)
    maximum = 0.0
    for fixed_four in combinations(range(dimension), 4):
        fixed_set = set(fixed_four)
        for singleton in range(dimension):
            energy = 0.0
            for element in range(dimension):
                if element in fixed_set:
                    continue
                support = tuple(sorted(fixed_four + (element,)))
                moment = (
                    quintic_weight(support, order)
                    * hadamard[singleton, support_xor(support)]
                )
                energy += moment * moment
            maximum = max(maximum, energy)
    return maximum


def exact_q4_base_chain() -> tuple[int, float, float]:
    """Check rank--Frobenius and an attaining fixed-cubic law."""

    order = 4
    dimension = order * order
    triples = tuple(combinations(range(dimension), 3))
    hadamard = sylvester(dimension)
    matrix = np.zeros((dimension, len(triples) * dimension))
    for triple_index, triple in enumerate(triples):
        weight = cubic_weight(triple, order)
        xor_triple = support_xor(triple)
        for b in range(dimension):
            column = triple_index * dimension + b
            matrix[:, column] = (
                weight * hadamard[xor_triple, b] * hadamard[b, :]
            )
    rank = int(np.linalg.matrix_rank(matrix))
    maximum_entry = float(np.max(np.abs(matrix)))

    vertical = (0, order, 2 * order)
    vertical_index = triples.index(vertical)
    block = matrix[
        :, vertical_index * dimension : (vertical_index + 1) * dimension
    ]
    weighted_block = block / dimension
    attained = float(np.linalg.svd(weighted_block, compute_uv=False).sum())
    return rank, maximum_entry, attained


def direct_q4_sparse_chain_checks() -> float:
    """Stress the exact indexed tensor under correlated diagonal laws."""

    rng = np.random.default_rng(SEED)
    order = 4
    dimension = order * order
    triples = tuple(combinations(range(dimension), 3))
    fours = tuple(combinations(range(dimension), 4))
    hadamard = sylvester(dimension)
    bound = column_cubic_quintic_coefficient(dimension)
    worst = 0.0

    all_rows = tuple(
        (singleton, fixed_four)
        for singleton in range(dimension)
        for fixed_four in fours
    )
    all_columns = tuple(
        (triple, b, element)
        for triple in triples
        for b in range(dimension)
        for element in range(dimension)
    )
    for _ in range(12):
        row_indices = rng.choice(len(all_rows), size=28, replace=False)
        column_indices = rng.choice(
            len(all_columns), size=44, replace=False
        )
        rows = tuple(all_rows[index] for index in row_indices)
        columns = tuple(all_columns[index] for index in column_indices)
        tensor = np.zeros((len(rows), len(columns)))
        for row_index, (c, fixed_four) in enumerate(rows):
            fixed_set = set(fixed_four)
            for column_index, (triple, b, element) in enumerate(columns):
                if element in fixed_set:
                    continue
                five = tuple(sorted(fixed_four + (element,)))
                tensor[row_index, column_index] = (
                    cubic_weight(triple, order)
                    * hadamard[support_xor(triple), b]
                    * hadamard[b, c]
                    * quintic_weight(five, order)
                    * hadamard[c, support_xor(five)]
                )
        row_law = rng.dirichlet(np.ones(len(rows)))
        column_law = rng.dirichlet(np.ones(len(columns)))
        weighted = (
            np.sqrt(row_law)[:, None]
            * tensor
            * np.sqrt(column_law)[None, :]
        )
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > bound * (1 + 4e-12):
            raise AssertionError(("sparse exact target", nuclear, bound))
        worst = max(worst, nuclear / bound)
    return worst


def main() -> None:
    q4_energy = exact_q4_quintic_row_energy()
    expected_q4_energy = endpoint_quintic_singleton_slice_energies(4)[4]
    if not np.isclose(q4_energy, expected_q4_energy, atol=2e-14):
        raise AssertionError(("q4 quintic row energy", q4_energy))

    rank, maximum_entry, attained = exact_q4_base_chain()
    if rank != 16:
        raise AssertionError(("base rank", rank))
    if not np.isclose(maximum_entry, 1 / 16, atol=2e-15):
        raise AssertionError(("base maximum entry", maximum_entry))
    if not np.isclose(attained, 1 / 4, atol=2e-15):
        raise AssertionError(("base rank-Frobenius saturation", attained))

    orbit = column_cubic_quintic_orbit_entries()
    if len(orbit) != 4:
        raise AssertionError(("column-cubic quintic orbit", orbit))
    coefficient = column_cubic_quintic_coefficient()
    if not np.isclose(coefficient, 0.031188905122404905, atol=2e-15):
        raise AssertionError(("row-energy coefficient", coefficient))
    result = column_cubic_quintic_contraction()
    if not result.generic_disjointness_coefficient > result.acceptance_gate:
        raise AssertionError(("generic factor should fail", result))
    if not result.coefficient < result.provisional_coefficient:
        raise AssertionError(("provisional improvement", result))
    if not result.coefficient < result.acceptance_gate:
        raise AssertionError(("acceptance gate", result))
    if not result.threshold_slack > 0.000967:
        raise AssertionError(("five-theorem ledger", result))
    if result.next_unresolved_entries[0] != (
        (1, 1, 3, 5),
        (0, 1, 1, 3),
    ):
        raise AssertionError(("reranked next orbit", result))
    if not np.isclose(
        result.next_admissible_coefficient,
        0.0570749885141564,
        atol=4e-14,
    ):
        raise AssertionError(("next admissible coefficient", result))

    worst = direct_q4_sparse_chain_checks()
    print(
        "column-cubic quintic-row contraction passed: "
        f"q4_row_energy={q4_energy:.12g},"
        f"q4_base_rank={rank},"
        f"q4_base_attained={attained:.12g},"
        f"q4_sparse_ratio={worst:.12g},"
        f"coefficient={result.coefficient:.12g},"
        f"ledger_total={result.optimized_total:.12g},"
        f"threshold_slack={result.threshold_slack:.12g},"
        f"next_admissible={result.next_admissible_coefficient:.12g}"
    )


if __name__ == "__main__":
    main()
