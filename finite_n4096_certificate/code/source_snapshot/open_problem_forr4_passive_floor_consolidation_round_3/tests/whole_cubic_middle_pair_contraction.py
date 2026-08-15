#!/usr/bin/env python3
"""Regression for the ninth balanced whole-cubic/middle-pair theorem."""

from __future__ import annotations

from itertools import combinations
from math import comb
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
from opposite_endpoint_orbit_scan import sylvester  # noqa: E402
from whole_cubic_middle_pair_contraction import (  # noqa: E402
    record_one_pair_incidence,
    record_sector_bounds,
    record_three_pair_counts,
    whole_cubic_middle_pair_coefficient,
    whole_cubic_middle_pair_contraction,
    whole_cubic_middle_pair_orbit_entries,
)


SEED = 2026071509


def row_multiplicities(
    order: int,
    support: tuple[int, ...],
) -> tuple[int, ...]:
    """Return the nonzero physical-row multiplicities."""

    counts: dict[int, int] = {}
    for cell in support:
        row = cell // order
        counts[row] = counts.get(row, 0) + 1
    return tuple(sorted(counts.values(), reverse=True))


def exact_pair_record_counts(
    order: int,
    fixed_pair: tuple[int, int],
) -> tuple[int, int, int]:
    """Enumerate record-one, 3+1+1, and 2+1+1+1 completions."""

    dimension = order * order
    fixed = set(fixed_pair)
    available = tuple(cell for cell in range(dimension) if cell not in fixed)
    record_one = 0
    record_three_no_even = 0
    record_three_one_even = 0
    for extra in combinations(available, 3):
        support = tuple(sorted(fixed_pair + extra))
        record = parity_record_size(order, support, axis=0)
        if record == 1:
            record_one += 1
        elif record == 3:
            pattern = row_multiplicities(order, support)
            if pattern == (3, 1, 1):
                record_three_no_even += 1
            elif pattern == (2, 1, 1, 1):
                record_three_one_even += 1
            else:
                raise AssertionError(("record-three row pattern", support, pattern))
    return record_one, record_three_no_even, record_three_one_even


def count_formula_checks() -> None:
    """Verify both fixed-pair geometries through q=8."""

    for order in range(4, 9):
        same = exact_pair_record_counts(order, (0, 1))
        distinct = exact_pair_record_counts(order, (0, order))
        counts = record_three_pair_counts(order)
        if same[0] != record_one_pair_incidence(order):
            raise AssertionError(("record-one same-row count", order, same))
        if not same[0] > distinct[0]:
            raise AssertionError(
                ("record-one maximizing geometry", order, same, distinct)
            )
        if same[1:] != (
            counts.same_row_no_even,
            counts.same_row_one_even,
        ):
            raise AssertionError(("record-three same-row counts", order, same, counts))
        if distinct[1:] != (
            counts.distinct_rows_no_even,
            counts.distinct_rows_one_even,
        ):
            raise AssertionError(
                ("record-three distinct-row counts", order, distinct, counts)
            )


def exact_q4_slice_table() -> tuple[float, float, float, float]:
    """Construct every q=4 M_35 row through every fixed pair."""

    order = 4
    dimension = order * order
    moments = exact_link_moments(order)
    pairs = tuple(combinations(range(dimension), 2))
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    incidence = sparse.lil_matrix((len(moments.supports_five), len(pairs)))
    for quintic_index, quintic in enumerate(moments.supports_five):
        for pair in combinations(quintic, 2):
            incidence[quintic_index, pair_index[pair]] = 1
    slice_table = np.square(moments.moment_35) @ incidence.tocsr()

    endpoint_energy = np.square(moments.moment_13).sum(axis=0)
    record_one_rows = np.asarray(
        [
            parity_record_size(order, cubic, axis=1) == 1 and endpoint_energy[index] > 0
            for index, cubic in enumerate(moments.supports_three)
        ]
    )
    record_three_rows = np.asarray(
        [
            parity_record_size(order, cubic, axis=1) == 3 and endpoint_energy[index] > 0
            for index, cubic in enumerate(moments.supports_three)
        ]
    )
    record_one_maximum = float(slice_table[record_one_rows].max())
    record_three_maximum = float(slice_table[record_three_rows].max())
    weighted_one = float(
        (endpoint_energy[:, None] * slice_table)[record_one_rows].max()
    )
    weighted_three = float(
        (endpoint_energy[:, None] * slice_table)[record_three_rows].max()
    )
    return record_one_maximum, record_three_maximum, weighted_one, weighted_three


def exact_representative_slice(
    order: int,
    cubic: tuple[int, int, int],
    fixed_pair: tuple[int, int],
) -> tuple[float, float, float, float]:
    """Enumerate one fixed cubic/pair M_35 squared slice by row pattern."""

    dimension = order * order
    fixed = set(fixed_pair)
    available = tuple(cell for cell in range(dimension) if cell not in fixed)
    no_even = 0.0
    one_even = 0.0
    no_even_maximum = 0.0
    one_even_maximum = 0.0
    for extra in combinations(available, 3):
        quintic = tuple(sorted(fixed_pair + extra))
        moment = combined_link_moment(order, cubic, quintic)
        pattern = row_multiplicities(order, quintic)
        if pattern == (3, 1, 1):
            no_even += moment**2
            no_even_maximum = max(no_even_maximum, abs(moment))
        elif pattern == (2, 1, 1, 1):
            one_even += moment**2
            one_even_maximum = max(one_even_maximum, abs(moment))
        elif moment != 0:
            raise AssertionError(("nonzero incompatible row pattern", quintic, moment))
    return no_even, one_even, no_even_maximum, one_even_maximum


def q8_slice_checks() -> tuple[float, float, float]:
    """Check leading representatives and the two record-three entry bounds."""

    q = 8
    record_one_cubic = (0, q - 1, 2 * q)
    record_one_pair = (0, q - 1)
    record_one_slice = sum(
        combined_link_moment(
            q,
            record_one_cubic,
            tuple(sorted(record_one_pair + extra)),
        )
        ** 2
        for extra in combinations(
            tuple(cell for cell in range(q * q) if cell not in set(record_one_pair)),
            3,
        )
    )
    if not np.isclose(record_one_slice, 1.170153061224345, atol=3e-12):
        raise AssertionError(("q8 record-one slice", record_one_slice))

    record_three_cubic = (0, 1, 2)
    same = exact_representative_slice(q, record_three_cubic, (0, q - 1))
    distinct = exact_representative_slice(q, record_three_cubic, (0, 2 * q))
    if not np.isclose(sum(same[:2]), 0.46666666666615036, atol=3e-12):
        raise AssertionError(("q8 same-row record-three slice", same))
    if not np.isclose(sum(distinct[:2]), 0.25000000000001527, atol=3e-12):
        raise AssertionError(("q8 distinct-row record-three slice", distinct))
    no_even_bound = 1 / comb(q, 3)
    one_even_bound = 3 / ((q - 3) * 56)
    for name, values in (("same", same), ("distinct", distinct)):
        if values[2] > no_even_bound * (1 + 2e-12):
            raise AssertionError(("q8 no-even entry maximum", name, values))
        if values[3] > one_even_bound * (1 + 2e-12):
            raise AssertionError(("q8 one-even entry maximum", name, values))
    return record_one_slice, sum(same[:2]), sum(distinct[:2])


def residual_walsh_check() -> float:
    """Stress the residual H M_13 matrix under correlated q=4 laws."""

    rng = np.random.default_rng(SEED)
    q = 4
    dimension = q * q
    moments = exact_link_moments(q)
    hadamard = sylvester(dimension)
    matrix = np.zeros((len(moments.supports_three), dimension * dimension))
    for cubic_index in range(len(moments.supports_three)):
        for a in range(dimension):
            for b in range(dimension):
                matrix[cubic_index, a * dimension + b] = (
                    hadamard[a, b] * moments.moment_13[b, cubic_index]
                )
    worst = 0.0
    for _ in range(24):
        row = rng.dirichlet(np.ones(len(matrix)))
        column = rng.dirichlet(np.ones(matrix.shape[1]))
        weighted = np.sqrt(row)[:, None] * matrix * np.sqrt(column)[None, :]
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > 1 / q * (1 + 3e-12):
            raise AssertionError(("residual Walsh compression", nuclear))
        worst = max(worst, q * nuclear)
    return worst


def direct_sparse_tensor_checks(coefficient: float) -> float:
    """Stress exact q=4 target submatrices under correlated laws."""

    rng = np.random.default_rng(SEED + 1)
    q = 4
    dimension = q * q
    triples = tuple(combinations(range(dimension), 3))
    pairs = tuple(combinations(range(dimension), 2))
    hadamard = sylvester(dimension)
    worst = 0.0
    for _ in range(12):
        rows = tuple(
            (
                triples[int(rng.integers(len(triples)))],
                pairs[int(rng.integers(len(pairs)))],
            )
            for _ in range(28)
        )
        columns = tuple(
            (
                int(rng.integers(dimension)),
                int(rng.integers(dimension)),
                triples[int(rng.integers(len(triples)))],
            )
            for _ in range(42)
        )
        tensor = np.zeros((len(rows), len(columns)))
        for row_index, (cubic, fixed_pair) in enumerate(rows):
            fixed = set(fixed_pair)
            for column_index, (first, second, extra) in enumerate(columns):
                if fixed.intersection(extra):
                    continue
                quintic = tuple(sorted(fixed_pair + extra))
                tensor[row_index, column_index] = (
                    hadamard[first, second]
                    * record_one_link_moment(q, (second,), cubic)
                    * combined_link_moment(q, cubic, quintic)
                )
        row_law = rng.dirichlet(np.ones(len(rows)))
        column_law = rng.dirichlet(np.ones(len(columns)))
        weighted = np.sqrt(row_law)[:, None] * tensor * np.sqrt(column_law)[None, :]
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > coefficient * (1 + 5e-12):
            raise AssertionError(("sparse target tensor", nuclear, coefficient))
        worst = max(worst, nuclear / coefficient)
    return worst


def main() -> None:
    count_formula_checks()

    q4_one, q4_three, q4_weighted_one, q4_weighted_three = exact_q4_slice_table()
    if not np.isclose(q4_one, 2.75000000000002, atol=4e-14):
        raise AssertionError(("q4 record-one slice", q4_one))
    if not np.isclose(q4_three, 2.5000000000000098, atol=4e-14):
        raise AssertionError(("q4 record-three slice", q4_three))
    if not np.isclose(q4_weighted_one, 0.30555555555555786, atol=4e-14):
        raise AssertionError(("q4 weighted record one", q4_weighted_one))
    if not np.isclose(q4_weighted_three, 2.5000000000000098, atol=4e-14):
        raise AssertionError(("q4 weighted record three", q4_weighted_three))

    q8_one, q8_same, q8_distinct = q8_slice_checks()
    residual_ratio = residual_walsh_check()

    bounds = record_sector_bounds(32)
    if bounds[0] != 15811580:
        raise AssertionError(("q32 record-one incidence", bounds[0]))
    if not np.isclose(bounds[2], 20.6379829857498, atol=3e-13):
        raise AssertionError(("q32 record-one slice bound", bounds[2]))
    if not np.isclose(bounds[3], 0.004579541014076183, atol=3e-16):
        raise AssertionError(("q32 record-one coefficient", bounds[3]))
    if not np.isclose(bounds[5], 0.6447163515016685, atol=3e-15):
        raise AssertionError(("q32 record-three same-row bound", bounds[5]))
    if not np.isclose(bounds[6], 0.060358731924360406, atol=3e-16):
        raise AssertionError(("q32 record-three distinct-row bound", bounds[6]))
    coefficient = whole_cubic_middle_pair_coefficient()
    if not np.isclose(coefficient, 0.025091947154681882, atol=3e-16):
        raise AssertionError(("q32 ninth coefficient", coefficient))

    result = whole_cubic_middle_pair_contraction()
    if not result.coefficient < result.provisional_coefficient:
        raise AssertionError(("provisional improvement", result))
    if not result.coefficient < result.acceptance_gate:
        raise AssertionError(("ninth-orbit gate", result))
    if len(whole_cubic_middle_pair_orbit_entries()) != 4:
        raise AssertionError(
            ("ninth target orbit", whole_cubic_middle_pair_orbit_entries())
        )
    if not np.isclose(result.optimized_total, 0.3326862124340385, atol=4e-11):
        raise AssertionError(("nine-theorem ledger", result))
    if not result.threshold_slack > 0.000647:
        raise AssertionError(("nine-theorem slack", result))
    if result.next_unresolved_entries[0] != (
        (1, 3, 5, 1),
        (0, 2, 3, 0),
    ):
        raise AssertionError(("reranked next orbit", result))
    if not np.isclose(
        result.next_admissible_coefficient,
        0.0529177166679669,
        atol=4e-14,
    ):
        raise AssertionError(("next coefficient gate", result))

    sparse_ratio = direct_sparse_tensor_checks(record_sector_bounds(4)[7])
    print(
        "whole-cubic middle-pair contraction passed: "
        f"q4_record_one={q4_one:.12g},"
        f"q4_record_three={q4_three:.12g},"
        f"q8_record_one={q8_one:.12g},"
        f"q8_record_three_same={q8_same:.12g},"
        f"q8_record_three_distinct={q8_distinct:.12g},"
        f"residual_ratio={residual_ratio:.12g},"
        f"coefficient={result.coefficient:.12g},"
        f"ledger_total={result.optimized_total:.12g},"
        f"threshold_slack={result.threshold_slack:.12g},"
        f"next_admissible={result.next_admissible_coefficient:.12g},"
        f"sparse_ratio={sparse_ratio:.12g}"
    )


if __name__ == "__main__":
    main()
