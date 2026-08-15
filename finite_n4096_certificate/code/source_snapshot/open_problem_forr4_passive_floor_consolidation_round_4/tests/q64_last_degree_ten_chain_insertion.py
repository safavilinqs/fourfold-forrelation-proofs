#!/usr/bin/env python3
"""Regression for the final q64 degree-ten chain theorem."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_last_degree_ten_chain_insertion import (  # noqa: E402
    TARGET,
    artifact_text,
    diagnostic,
    inserted_coefficients,
    last_degree_ten_coefficient,
    last_degree_ten_entries,
    pre_last_degree_ten_quintic_entries,
    remaining_quintic_entries,
)
from q64_whole_cubic_decorated_row_insertion import (  # noqa: E402
    inserted_coefficients as previous_inserted_coefficients,
)
from adjacent_cubic_quintic_orbit_witness import (  # noqa: E402
    record_one_link_moment,
)
from opposite_endpoint_orbit_scan import (  # noqa: E402
    endpoint_moment,
    sylvester,
)


def sparse_exact_chain_stress(order: int) -> float:
    """Stress exact target submatrices under correlated diagonal laws."""

    rng = np.random.default_rng(2026071612 + order)
    dimension = order * order
    hadamard = sylvester(dimension)
    pairs = tuple(combinations(range(dimension), 2))
    triples = tuple(combinations(range(dimension), 3))
    coefficient = last_degree_ten_coefficient(order)
    worst = 0.0
    for _ in range(24):
        rows = tuple(
            (
                int(rng.integers(dimension)),
                int(rng.integers(dimension)),
                triples[int(rng.integers(len(triples)))],
            )
            for _ in range(24)
        )
        columns = tuple(
            (
                pairs[int(rng.integers(len(pairs)))],
                int(rng.integers(dimension)),
                pairs[int(rng.integers(len(pairs)))],
            )
            for _ in range(36)
        )
        matrix = np.zeros((len(rows), len(columns)))
        for row_index, (cell, first_singleton, fixed_triple) in enumerate(rows):
            triple_set = set(fixed_triple)
            for column_index, (fixed_pair, second_singleton, added_pair) in enumerate(
                columns
            ):
                if cell in fixed_pair or triple_set.intersection(added_pair):
                    continue
                cubic = tuple(sorted(fixed_pair + (cell,)))
                quintic = tuple(sorted(fixed_triple + added_pair))
                matrix[row_index, column_index] = (
                    record_one_link_moment(
                        order,
                        (first_singleton,),
                        cubic,
                    )
                    * hadamard[first_singleton, second_singleton]
                    * endpoint_moment(
                        quintic,
                        second_singleton,
                        order,
                        5,
                        False,
                    )
                )
        row_law = rng.dirichlet(np.ones(len(rows)))
        column_law = rng.dirichlet(np.ones(len(columns)))
        weighted = (
            np.sqrt(row_law)[:, None]
            * matrix
            * np.sqrt(column_law)[None, :]
        )
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > coefficient * (1 + 5e-12):
            raise AssertionError(
                ("final degree-ten sparse tensor", order, nuclear, coefficient)
            )
        worst = max(worst, nuclear / coefficient)
    return worst


def walsh_chain_residual_stress(order: int) -> float:
    """Stress the exact ``1/q`` residual after endpoint completion."""

    rng = np.random.default_rng(2026071712 + order)
    dimension = order * order
    hadamard = sylvester(dimension)
    worst = 0.0
    for _ in range(32):
        row_labels = rng.integers(dimension, size=26)
        column_labels = rng.integers(dimension, size=40)
        matrix = hadamard[np.ix_(row_labels, column_labels)] / order
        row_law = rng.dirichlet(np.ones(len(row_labels)))
        column_law = rng.dirichlet(np.ones(len(column_labels)))
        weighted = (
            np.sqrt(row_law)[:, None]
            * matrix
            * np.sqrt(column_law)[None, :]
        )
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > 1 / order + 5e-12:
            raise AssertionError(("final degree-ten Walsh residual", nuclear))
        worst = max(worst, nuclear)
    return worst


def main() -> None:
    result = diagnostic()
    discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.closed_entries,
        result.closed_orbits,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
        result.remaining_quintic_entries,
        result.remaining_extreme_entries,
        result.remaining_balanced_entries,
    )
    expected_discrete = (64, 4096, 16_384, 4, 1, 320, 324, 564, 104, 72, 32)
    if discrete != expected_discrete:
        raise AssertionError(("last degree-ten discrete result", discrete))

    observed = (
        result.cubic_distinctness_factor,
        result.quintic_completion_factor,
        result.walsh_chain_factor,
        result.coefficient,
        result.previous_routing.total,
        result.last_degree_ten_inserted.total,
        result.last_degree_ten_inserted.beta,
        result.last_degree_ten_inserted.perron_upper,
        result.last_degree_ten_inserted.promise_loss,
        result.last_degree_ten_inserted.margin_to_one_third,
        result.routing_margin_improvement,
        result.remaining_quintic_local_proxy.total,
        result.remaining_quintic_local_proxy.margin_to_one_third,
        result.proxy_reserve_after_declared_allowance,
    )
    expected = (
        2.413207370108021,
        2.414213562373095,
        0.015625,
        0.09103121815208583,
        0.32528397960750105,
        0.3247663269533303,
        0.7461529246322275,
        0.3077302268911228,
        0.017036100062207536,
        0.008567006380003017,
        0.0005176526541707527,
        0.3250376129305583,
        0.008295720402775042,
        0.007295720402775042,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("last degree-ten numeric result", observed))
    if not np.isclose(
        result.coefficient,
        result.cubic_distinctness_factor
        * result.quintic_completion_factor
        * result.walsh_chain_factor,
        rtol=1e-13,
    ):
        raise AssertionError("last degree-ten coefficient identity")

    entries = set(last_degree_ten_entries())
    pre_last = set(pre_last_degree_ten_quintic_entries())
    remaining = set(remaining_quintic_entries())
    if TARGET not in entries or len(entries) != 4:
        raise AssertionError(("last degree-ten orbit inventory", entries))
    if not entries.issubset(pre_last) or entries.intersection(remaining):
        raise AssertionError("last degree-ten partition overlap")
    if len(entries) + len(remaining) != len(pre_last):
        raise AssertionError("last degree-ten partition size")
    for entry in entries:
        if np.isclose(
            inserted_coefficients()[entry],
            previous_inserted_coefficients()[entry],
        ):
            raise AssertionError(("last degree-ten theorem not inserted", entry))

    sparse_worst = max(sparse_exact_chain_stress(order) for order in (4, 8))
    residual_worst = max(walsh_chain_residual_stress(order) for order in (4, 8))

    committed = (
        ROOT / "artifacts" / "q64_last_degree_ten_chain_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 final degree-ten artifact")

    print(
        "q64 final degree-ten chain insertion passed: "
        f"entries={result.closed_entries},"
        f"coefficient={result.coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.last_degree_ten_inserted.total:.12g},"
        f"margin={result.last_degree_ten_inserted.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        f"sparse_stress={sparse_worst:.12g},"
        f"residual_stress={residual_worst:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
