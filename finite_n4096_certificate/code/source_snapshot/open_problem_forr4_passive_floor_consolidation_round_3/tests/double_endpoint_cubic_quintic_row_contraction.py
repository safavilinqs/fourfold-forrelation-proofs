#!/usr/bin/env python3
"""Regression for the tenth balanced double-endpoint row theorem."""

from __future__ import annotations

from functools import reduce
from itertools import combinations
from operator import xor
from pathlib import Path
import sys

import numpy as np
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from adjacent_cubic_quintic_mixed_orbit_q4 import (  # noqa: E402
    combined_link_moment,
)
from adjacent_cubic_quintic_orbit_witness import (  # noqa: E402
    exact_link_moments,
    parity_record_size,
    record_one_link_moment,
)
from double_endpoint_cubic_quintic_row_contraction import (  # noqa: E402
    double_endpoint_cubic_quintic_coefficient,
    double_endpoint_cubic_quintic_contraction,
    double_endpoint_cubic_quintic_orbit_entries,
    scalar_row_sector_bounds,
)
from occupation_compatible_sector_optimization import (  # noqa: E402
    endpoint_quintic_singleton_slice_energies,
    endpoint_singleton_slice_energies,
    middle_quintic_incidence_bound,
)
from opposite_endpoint_orbit_scan import (  # noqa: E402
    cubic_weight,
    endpoint_moment,
    quintic_weight,
    sylvester,
)


SEED = 2026071510


def transpose_support(
    order: int,
    support: tuple[int, ...],
) -> tuple[int, ...]:
    """Swap the physical row and column labels of a support."""

    return tuple(sorted((cell % order) * order + cell // order for cell in support))


def exact_q4_scalar_row_table() -> tuple[
    float,
    float,
    float,
    tuple[int, int],
    tuple[int, int, int],
]:
    """Compute every q=4 scalar completion-row energy."""

    q = 4
    dimension = q * q
    moments = exact_link_moments(q)
    pairs = tuple(combinations(range(dimension), 2))
    triples = moments.supports_three
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    triple_index = {triple: index for index, triple in enumerate(triples)}

    cubic_weights = np.asarray(
        [cubic_weight(transpose_support(q, cubic), q) for cubic in triples]
    )
    quintic_weights = np.asarray(
        [quintic_weight(quintic, q) for quintic in moments.supports_five]
    )
    weighted = np.square(moments.moment_35)
    weighted *= np.square(cubic_weights)[:, None]
    weighted *= np.square(quintic_weights)[None, :]

    cubic_incidence = sparse.lil_matrix((len(pairs), len(triples)))
    for cubic_position, cubic in enumerate(triples):
        for pair in combinations(cubic, 2):
            cubic_incidence[pair_index[pair], cubic_position] = 1
    quintic_incidence = sparse.lil_matrix((len(moments.supports_five), len(triples)))
    for quintic_position, quintic in enumerate(moments.supports_five):
        for triple in combinations(quintic, 3):
            quintic_incidence[quintic_position, triple_index[triple]] = 1

    records = np.asarray([parity_record_size(q, cubic, axis=1) for cubic in triples])
    record_tables = {}
    for record in (1, 3):
        sector = weighted * (records == record)[:, None]
        record_tables[record] = (
            cubic_incidence.tocsr() @ sector @ quintic_incidence.tocsr()
        )
    total = record_tables[1] + record_tables[3]
    maximum_index = np.unravel_index(int(np.argmax(total)), total.shape)
    return (
        float(total[maximum_index]),
        float(record_tables[1].max()),
        float(record_tables[3].max()),
        pairs[maximum_index[0]],
        triples[maximum_index[1]],
    )


def exact_q8_representative_row() -> tuple[float, float]:
    """Evaluate the horizontal-pair/vertical-triple representative."""

    q = 8
    dimension = q * q
    fixed_pair = (0, 1)
    fixed_triple = (q - 1, 3 * q - 1, 4 * q - 1)
    pair_set = set(fixed_pair)
    triple_set = set(fixed_triple)
    energies = {1: 0.0, 3: 0.0}
    for cell in range(dimension):
        if cell in pair_set:
            continue
        cubic = tuple(sorted(fixed_pair + (cell,)))
        cubic_factor = cubic_weight(transpose_support(q, cubic), q) ** 2
        if cubic_factor == 0:
            continue
        record = parity_record_size(q, cubic, axis=1)
        for added in combinations(
            tuple(value for value in range(dimension) if value not in triple_set),
            2,
        ):
            quintic = tuple(sorted(fixed_triple + added))
            quintic_factor = quintic_weight(quintic, q) ** 2
            if quintic_factor == 0:
                continue
            moment = combined_link_moment(q, cubic, quintic)
            energies[record] += cubic_factor * moment**2 * quintic_factor
    return energies[1], energies[3]


def residual_walsh_tensor_check() -> float:
    """Stress the repeated, column-twisted H_N tensor H_N base."""

    rng = np.random.default_rng(SEED)
    q = 4
    dimension = q * q
    hadamard = sylvester(dimension)
    pairs = tuple(combinations(range(dimension), 2))
    triples = tuple(combinations(range(dimension), 3))
    worst = 0.0
    for _ in range(24):
        rows = tuple(
            (
                pairs[int(rng.integers(len(pairs)))],
                triples[int(rng.integers(len(triples)))],
            )
            for _ in range(30)
        )
        columns = tuple(
            (
                int(rng.integers(dimension)),
                int(rng.integers(dimension)),
                pairs[int(rng.integers(len(pairs)))],
                int(rng.integers(dimension)),
            )
            for _ in range(44)
        )
        matrix = np.zeros((len(rows), len(columns)))
        for row_index, (fixed_pair, fixed_triple) in enumerate(rows):
            pair_xor = reduce(xor, fixed_pair, 0)
            triple_xor = reduce(xor, fixed_triple, 0)
            for column_index, (first, cell, added_pair, final) in enumerate(columns):
                matrix[row_index, column_index] = (
                    hadamard[pair_xor ^ cell, first]
                    * hadamard[triple_xor ^ reduce(xor, added_pair, 0), final]
                )
        row_law = rng.dirichlet(np.ones(len(rows)))
        column_law = rng.dirichlet(np.ones(len(columns)))
        weighted = np.sqrt(row_law)[:, None] * matrix * np.sqrt(column_law)[None, :]
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > 1 + 5e-12:
            raise AssertionError(("residual Walsh tensor", nuclear))
        worst = max(worst, nuclear)
    return worst


def direct_sparse_tensor_checks(q4_coefficient: float) -> float:
    """Stress exact q=4 target submatrices under correlated laws."""

    rng = np.random.default_rng(SEED + 1)
    q = 4
    dimension = q * q
    pairs = tuple(combinations(range(dimension), 2))
    triples = tuple(combinations(range(dimension), 3))
    worst = 0.0
    for _ in range(12):
        rows = tuple(
            (
                pairs[int(rng.integers(len(pairs)))],
                triples[int(rng.integers(len(triples)))],
            )
            for _ in range(26)
        )
        columns = tuple(
            (
                int(rng.integers(dimension)),
                int(rng.integers(dimension)),
                pairs[int(rng.integers(len(pairs)))],
                int(rng.integers(dimension)),
            )
            for _ in range(40)
        )
        tensor = np.zeros((len(rows), len(columns)))
        for row_index, (fixed_pair, fixed_triple) in enumerate(rows):
            pair_set = set(fixed_pair)
            triple_set = set(fixed_triple)
            for column_index, (first, cell, added_pair, final) in enumerate(columns):
                if cell in pair_set or triple_set.intersection(added_pair):
                    continue
                cubic = tuple(sorted(fixed_pair + (cell,)))
                quintic = tuple(sorted(fixed_triple + added_pair))
                tensor[row_index, column_index] = (
                    record_one_link_moment(q, (first,), cubic)
                    * combined_link_moment(q, cubic, quintic)
                    * endpoint_moment(quintic, final, q, 5, False)
                )
        row_law = rng.dirichlet(np.ones(len(rows)))
        column_law = rng.dirichlet(np.ones(len(columns)))
        weighted = np.sqrt(row_law)[:, None] * tensor * np.sqrt(column_law)[None, :]
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > q4_coefficient * (1 + 5e-12):
            raise AssertionError(("sparse target tensor", nuclear, q4_coefficient))
        worst = max(worst, nuclear / q4_coefficient)
    return worst


def main() -> None:
    q4_total, q4_one, q4_three, maximizing_pair, maximizing_triple = (
        exact_q4_scalar_row_table()
    )
    if not np.isclose(q4_total, 1.5869341563786017, atol=4e-14):
        raise AssertionError(("q4 total scalar row", q4_total))
    if not np.isclose(q4_one, 0.03360768175582992, atol=4e-14):
        raise AssertionError(("q4 record-one scalar row", q4_one))
    if not np.isclose(q4_three, 1.5781893004115233, atol=4e-14):
        raise AssertionError(("q4 record-three scalar row", q4_three))
    if maximizing_pair[0] // 4 != maximizing_pair[1] // 4:
        raise AssertionError(("q4 maximizing pair geometry", maximizing_pair))
    if len({cell % 4 for cell in maximizing_triple}) != 1:
        raise AssertionError(("q4 maximizing row", maximizing_pair, maximizing_triple))

    q8_one, q8_three = exact_q8_representative_row()
    if not np.isclose(q8_one, 4.1968057527050735e-06, atol=3e-16):
        raise AssertionError(("q8 record-one row", q8_one))
    if not np.isclose(q8_three, 0.005373186342425646, atol=3e-15):
        raise AssertionError(("q8 record-three row", q8_three))

    bounds = scalar_row_sector_bounds(32)
    expected_cubic = endpoint_singleton_slice_energies(32)[2]
    expected_quintic = endpoint_quintic_singleton_slice_energies(32)[3]
    if not np.isclose(bounds[0], expected_cubic, atol=2e-16):
        raise AssertionError(("q32 cubic pair slice", bounds))
    if not np.isclose(bounds[1], expected_quintic, atol=2e-15):
        raise AssertionError(("q32 quintic triple slice", bounds))
    if bounds[2] != 62 or bounds[3] != middle_quintic_incidence_bound(32, 3):
        raise AssertionError(("q32 record-one incidences", bounds))
    if bounds[3] != 3780:
        raise AssertionError(("q32 quintic triple incidence", bounds[3]))
    if not np.isclose(bounds[5], 6069 / 19066240, atol=2e-18):
        raise AssertionError(("q32 record-one bound", bounds[5]))
    if not np.isclose(bounds[7], 2151513 / 1182106880, atol=2e-18):
        raise AssertionError(("q32 record-three bound", bounds[7]))

    coefficient = double_endpoint_cubic_quintic_coefficient()
    if not np.isclose(coefficient, 0.04624259624455661, atol=3e-16):
        raise AssertionError(("q32 tenth coefficient", coefficient))
    result = double_endpoint_cubic_quintic_contraction()
    if not result.coefficient > result.provisional_coefficient:
        raise AssertionError(("expected provisional regression", result))
    if not result.coefficient < result.acceptance_gate:
        raise AssertionError(("tenth-orbit gate", result))
    if len(double_endpoint_cubic_quintic_orbit_entries()) != 4:
        raise AssertionError(
            ("tenth target orbit", double_endpoint_cubic_quintic_orbit_entries())
        )
    if not np.isclose(result.optimized_total, 0.333132605485488, atol=4e-11):
        raise AssertionError(("ten-theorem ledger", result))
    if not result.threshold_slack > 0.000200:
        raise AssertionError(("ten-theorem slack", result))
    if result.next_unresolved_entries[0] != (
        (1, 1, 5, 3),
        (0, 1, 3, 1),
    ):
        raise AssertionError(("reranked next orbit", result))
    if not np.isclose(
        result.next_admissible_coefficient, 0.0379251204234081, atol=4e-14
    ):
        raise AssertionError(("next coefficient gate", result))

    residual_worst = residual_walsh_tensor_check()
    sparse_worst = direct_sparse_tensor_checks(np.sqrt(q4_total))
    print(
        "double-endpoint cubic-quintic row contraction passed: "
        f"q4_row={q4_total:.12g},"
        f"q8_row={q8_one + q8_three:.12g},"
        f"q32_row={result.row_energy_bound:.12g},"
        f"coefficient={result.coefficient:.12g},"
        f"ledger_total={result.optimized_total:.12g},"
        f"threshold_slack={result.threshold_slack:.12g},"
        f"residual_worst={residual_worst:.12g},"
        f"sparse_ratio={sparse_worst:.12g}"
    )


if __name__ == "__main__":
    main()
