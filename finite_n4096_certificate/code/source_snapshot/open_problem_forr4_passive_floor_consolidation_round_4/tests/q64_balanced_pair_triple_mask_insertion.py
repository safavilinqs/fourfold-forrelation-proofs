#!/usr/bin/env python3
"""Regression for the q64 balanced pair--triple mask insertion."""

from __future__ import annotations

from itertools import combinations
from math import sqrt
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_balanced_pair_triple_mask_insertion import (  # noqa: E402
    artifact_text,
    balanced_pair_triple_entries,
    diagnostic,
    pair_triple_disjointness_factor,
    remaining_quintic_entries,
)
from q64_degree_ten_completion_row_insertion import orbit  # noqa: E402
from q64_dual_endpoint_schur_insertion import (  # noqa: E402
    has_favorable_cubic_singleton,
    has_favorable_quintic_singleton,
)
from q64_post_universal_quintic_gate import (  # noqa: E402
    quintic_split_depth,
)


def mask_factorization_check(dimension: int) -> None:
    pairs = tuple(combinations(range(dimension), 2))
    triples = tuple(combinations(range(dimension), 3))
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    x0 = 1.0
    x1 = sqrt(3 / 2)
    x2 = sqrt(3)
    factor = pair_triple_disjointness_factor()
    maximum_row_norm = 0.0
    maximum_column_norm = 0.0
    for pair in pairs:
        row = np.zeros(1 + dimension + len(pairs))
        row[0] = sqrt(x0)
        row[1 + np.asarray(pair)] = sqrt(x1)
        row[1 + dimension + pair_index[pair]] = sqrt(x2)
        maximum_row_norm = max(maximum_row_norm, np.linalg.norm(row))
        for triple in triples:
            column = np.zeros_like(row)
            column[0] = 1 / sqrt(x0)
            column[1 + np.asarray(triple)] = -1 / sqrt(x1)
            for contained in combinations(triple, 2):
                column[1 + dimension + pair_index[contained]] = 1 / sqrt(x2)
            expected = float(set(pair).isdisjoint(triple))
            if not np.isclose(np.dot(row, column), expected, atol=1e-14):
                raise AssertionError(("pair--triple mask", pair, triple))
            maximum_column_norm = max(
                maximum_column_norm, np.linalg.norm(column)
            )
    if not np.isclose(
        maximum_row_norm * maximum_column_norm,
        factor,
        rtol=1e-13,
        atol=1e-13,
    ):
        raise AssertionError(("pair--triple factor", dimension))


def main() -> None:
    result = diagnostic()
    entries = balanced_pair_triple_entries()
    if len(entries) != 8:
        raise AssertionError(("balanced pair--triple entries", len(entries)))
    if len({frozenset(orbit(entry)) for entry in entries}) != 2:
        raise AssertionError("balanced pair--triple orbit count")
    for entry in entries:
        if quintic_split_depth(entry) != 2:
            raise AssertionError(("nonbalanced entry", entry))
        if not has_favorable_cubic_singleton(entry):
            raise AssertionError(("missing cubic endpoint", entry))
        if has_favorable_quintic_singleton(entry):
            raise AssertionError(("unexpected dual endpoint", entry))
    mask_factorization_check(5)
    mask_factorization_check(7)
    discrete = (
        result.closed_entries,
        result.closed_orbits,
        result.extreme_entries,
        result.balanced_entries,
        result.quintic_favorable_entries,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
        result.remaining_quintic_entries,
        result.remaining_extreme_entries,
        result.remaining_balanced_entries,
        len(remaining_quintic_entries()),
    )
    if discrete != (8, 2, 0, 8, 0, 340, 348, 540, 80, 56, 24, 80):
        raise AssertionError(("balanced pair--triple discrete", discrete))
    observed = (
        result.cubic_fixed_pair_energy,
        result.cubic_endpoint_factor,
        result.pair_triple_mask_factor,
        result.coefficient,
        result.previous_routing.total,
        result.inserted_routing.total,
        result.inserted_routing.margin_to_one_third,
        result.routing_margin_spent,
        result.remaining_quintic_local_proxy.total,
        result.proxy_reserve_after_declared_allowance,
    )
    expected = (
        0.015384734623015874,
        0.12403521525357087,
        5.1815405503520555,
        0.6426934975083619,
        0.3247663269533303,
        0.32734300501757363,
        0.005990328315759685,
        0.0025766780642431653,
        0.3274960726492211,
        0.004837260684112222,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("balanced pair--triple numeric", observed))
    committed = (
        ROOT / "artifacts" / "q64_balanced_pair_triple_mask_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale balanced pair--triple artifact")
    print(
        "q64 balanced pair-triple mask insertion passed: "
        f"entries={result.closed_entries},"
        f"mask_factor={result.pair_triple_mask_factor:.12g},"
        f"coefficient={result.coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.inserted_routing.total:.12g},"
        f"margin={result.inserted_routing.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
