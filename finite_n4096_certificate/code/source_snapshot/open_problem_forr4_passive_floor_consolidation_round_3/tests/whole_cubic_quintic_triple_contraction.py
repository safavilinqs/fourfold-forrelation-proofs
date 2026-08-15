#!/usr/bin/env python3
"""Regression for the seventh balanced endpoint-slice contraction."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from occupation_compatible_sector_optimization import (  # noqa: E402
    endpoint_quintic_singleton_slice_energies,
)
from opposite_endpoint_orbit_scan import (  # noqa: E402
    cubic_weight,
    quintic_weight,
    support_xor,
    sylvester,
)
from whole_cubic_quintic_triple_contraction import (  # noqa: E402
    quintic_three_slice_energy,
    whole_cubic_quintic_triple_coefficient,
    whole_cubic_quintic_triple_contraction,
    whole_cubic_quintic_triple_orbit_entries,
)


SEED = 2026071507


def exact_quintic_three_slice(
    order: int,
    *,
    exhaustive: bool,
) -> tuple[float, tuple[int, int, int]]:
    """Evaluate exact M_15 fixed-triple rows from endpoint moments."""

    dimension = order * order
    hadamard = sylvester(dimension)
    if exhaustive:
        fixed_triples = tuple(combinations(range(dimension), 3))
    else:
        fixed_triples = ((0, order, 2 * order),)
    maximum = 0.0
    maximizing = fixed_triples[0]
    for fixed in fixed_triples:
        fixed_set = set(fixed)
        available = tuple(cell for cell in range(dimension) if cell not in fixed_set)
        for singleton in range(dimension):
            energy = 0.0
            for pair in combinations(available, 2):
                support = tuple(sorted(fixed + pair))
                moment = (
                    quintic_weight(support, order)
                    * hadamard[singleton, support_xor(support)]
                )
                energy += moment * moment
            if energy > maximum:
                maximum = energy
                maximizing = fixed
    return maximum, maximizing


def cubic_xor_compression_checks() -> float:
    """Check exact column compression and saturation of its unit bound."""

    rng = np.random.default_rng(SEED)
    order = 4
    dimension = order * order
    triples = tuple(combinations(range(dimension), 3))
    hadamard = sylvester(dimension)
    weights = np.asarray([cubic_weight(triple, order) for triple in triples])
    xors = np.asarray([support_xor(triple) for triple in triples])
    matrix = np.asarray(
        [
            [
                weights[index] * hadamard[xors[index], singleton]
                for index in range(len(triples))
            ]
            for singleton in range(dimension)
        ]
    )
    worst = 0.0
    for _ in range(16):
        row_law = rng.dirichlet(np.ones(dimension))
        column_law = rng.dirichlet(np.ones(len(triples)))
        weighted = np.sqrt(row_law)[:, None] * matrix * np.sqrt(column_law)[None, :]
        direct = float(np.linalg.svd(weighted, compute_uv=False).sum())
        xor_mass = np.bincount(
            xors,
            weights=column_law * weights * weights,
            minlength=dimension,
        )
        compressed = np.sqrt(row_law)[:, None] * hadamard.T * np.sqrt(xor_mass)[None, :]
        collapsed = float(np.linalg.svd(compressed, compute_uv=False).sum())
        if not np.isclose(direct, collapsed, atol=2e-13):
            raise AssertionError(("cubic xor compression", direct, collapsed))
        if direct > 1 + 2e-13:
            raise AssertionError(("cubic weighted trace norm", direct))
        worst = max(worst, direct)

    representatives = []
    for target_xor in range(dimension):
        candidates = [
            index
            for index, (xor_value, weight) in enumerate(zip(xors, weights, strict=True))
            if xor_value == target_xor and np.isclose(abs(weight), 1.0)
        ]
        if not candidates:
            raise AssertionError(("missing unit cubic xor", target_xor))
        representatives.append(candidates[0])
    row_law = np.full(dimension, 1 / dimension)
    column_law = np.zeros(len(triples))
    column_law[representatives] = 1 / dimension
    weighted = np.sqrt(row_law)[:, None] * matrix * np.sqrt(column_law)[None, :]
    attained = float(np.linalg.svd(weighted, compute_uv=False).sum())
    if not np.isclose(attained, 1.0, atol=2e-13):
        raise AssertionError(("cubic compression saturation", attained))
    return worst


def direct_sparse_tensor_checks() -> float:
    """Stress exact q=4 target submatrices under correlated laws."""

    rng = np.random.default_rng(SEED)
    order = 4
    dimension = order * order
    triples = tuple(combinations(range(dimension), 3))
    pairs = tuple(combinations(range(dimension), 2))
    hadamard = sylvester(dimension)
    bound = whole_cubic_quintic_triple_coefficient(dimension)
    rows = tuple(
        (b, c, fixed)
        for b in range(dimension)
        for c in range(dimension)
        for fixed in triples
    )
    columns = tuple((cubic, pair) for cubic in triples for pair in pairs)
    worst = 0.0
    for _ in range(12):
        selected_rows = tuple(
            rows[index] for index in rng.choice(len(rows), size=26, replace=False)
        )
        selected_columns = tuple(
            columns[index] for index in rng.choice(len(columns), size=40, replace=False)
        )
        tensor = np.zeros((len(selected_rows), len(selected_columns)))
        for row_index, (b, c, fixed) in enumerate(selected_rows):
            fixed_set = set(fixed)
            for column_index, (cubic, pair) in enumerate(selected_columns):
                if fixed_set.intersection(pair):
                    continue
                quintic = tuple(sorted(fixed + pair))
                tensor[row_index, column_index] = (
                    cubic_weight(cubic, order)
                    * hadamard[support_xor(cubic), b]
                    * hadamard[b, c]
                    * quintic_weight(quintic, order)
                    * hadamard[c, support_xor(quintic)]
                )
        row_law = rng.dirichlet(np.ones(len(selected_rows)))
        column_law = rng.dirichlet(np.ones(len(selected_columns)))
        weighted = np.sqrt(row_law)[:, None] * tensor * np.sqrt(column_law)[None, :]
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > bound * (1 + 5e-12):
            raise AssertionError(("sparse target tensor", nuclear, bound))
        worst = max(worst, nuclear / bound)
    return worst


def main() -> None:
    q4_energy, q4_maximizer = exact_quintic_three_slice(
        4,
        exhaustive=True,
    )
    if not np.isclose(q4_energy, 7 / 8, atol=2e-14):
        raise AssertionError(("q4 quintic triple slice", q4_energy))
    if len({cell % 4 for cell in q4_maximizer}) != 1:
        raise AssertionError(("q4 maximizing triple", q4_maximizer))
    q8_energy, _ = exact_quintic_three_slice(8, exhaustive=False)
    if not np.isclose(q8_energy, 261 / 224, atol=3e-14):
        raise AssertionError(("q8 vertical triple slice", q8_energy))

    for order in (4, 8, 16, 32):
        expected = endpoint_quintic_singleton_slice_energies(order)[3]
        observed = quintic_three_slice_energy(order)
        if not np.isclose(observed, expected, atol=3e-13):
            raise AssertionError(("closed quintic slice", order, observed))

    cubic_random_worst = cubic_xor_compression_checks()
    orbit = whole_cubic_quintic_triple_orbit_entries()
    if len(orbit) != 4:
        raise AssertionError(("seventh target orbit", orbit))
    coefficient = whole_cubic_quintic_triple_coefficient()
    if not np.isclose(coefficient, 0.037095279315720764, atol=3e-15):
        raise AssertionError(("q32 coefficient", coefficient))

    result = whole_cubic_quintic_triple_contraction()
    if not result.coefficient > result.provisional_coefficient:
        raise AssertionError(("safe coefficient should exceed provisional", result))
    if not result.coefficient < result.acceptance_gate:
        raise AssertionError(("seventh-orbit gate", result))
    if not np.isclose(
        result.optimized_total,
        0.33296436358896037,
        atol=4e-11,
    ):
        raise AssertionError(("seven-theorem ledger", result))
    if not result.threshold_slack > 0.000368:
        raise AssertionError(("seven-theorem slack", result))
    if result.next_unresolved_entries[0] != (
        (1, 3, 5, 1),
        (0, 2, 2, 1),
    ):
        raise AssertionError(("reranked next orbit", result))
    if not np.isclose(
        result.next_admissible_coefficient,
        0.04262693092906616,
        atol=4e-14,
    ):
        raise AssertionError(("next coefficient gate", result))

    sparse_worst = direct_sparse_tensor_checks()
    print(
        "whole-cubic quintic-triple contraction passed: "
        f"q4_slice={q4_energy:.12g},"
        f"q8_slice={q8_energy:.12g},"
        f"cubic_random_worst={cubic_random_worst:.12g},"
        f"coefficient={result.coefficient:.12g},"
        f"ledger_total={result.optimized_total:.12g},"
        f"threshold_slack={result.threshold_slack:.12g},"
        f"next_admissible={result.next_admissible_coefficient:.12g},"
        f"sparse_ratio={sparse_worst:.12g}"
    )


if __name__ == "__main__":
    main()
