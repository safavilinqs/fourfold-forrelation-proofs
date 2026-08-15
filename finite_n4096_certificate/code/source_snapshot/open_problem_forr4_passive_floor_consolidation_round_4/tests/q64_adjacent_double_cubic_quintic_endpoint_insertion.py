#!/usr/bin/env python3
"""Regression for the q64 adjacent-double-cubic endpoint insertion."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_adjacent_double_cubic_quintic_endpoint_insertion import (  # noqa: E402
    adjacent_double_cubic_entries,
    artifact_text,
    coefficient_map,
    diagnostic,
    endpoint_middle_split_pair,
    remaining_quintic_entries,
    TARGETS,
    target_coefficient,
    uniform_double_cubic_coefficient,
)
from q64_balanced_pair_triple_mask_insertion import (  # noqa: E402
    remaining_quintic_entries as pre_insertion_entries,
)
from q64_dual_endpoint_schur_insertion import (  # noqa: E402
    has_favorable_quintic_singleton,
)
from q64_post_universal_quintic_gate import (  # noqa: E402
    quintic_split_depth,
)


def main() -> None:
    result = diagnostic()
    entries = adjacent_double_cubic_entries()
    if len(entries) != 32 or not set(entries).issubset(pre_insertion_entries()):
        raise AssertionError(("adjacent double-cubic entries", entries))
    pairs = Counter(endpoint_middle_split_pair(entry) for entry in entries)
    expected_pairs = Counter(
        {
            pair: 4
            for pair in (
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 3),
                (2, 0),
                (2, 3),
                (3, 1),
                (3, 2),
            )
        }
    )
    if pairs != expected_pairs:
        raise AssertionError(("double-cubic split pairs", pairs))
    if sum(quintic_split_depth(entry) == 1 for entry in entries) != 24:
        raise AssertionError("adjacent double-cubic depth partition")
    if sum(has_favorable_quintic_singleton(entry) for entry in entries) != 12:
        raise AssertionError("favorable/generic quintic partition")
    q32 = tuple(
        uniform_double_cubic_coefficient(*pair, order=32)
        for pair in ((0, 1), (3, 2), (2, 2))
    )
    expected_q32 = (
        9.06838944474454e-06,
        0.0002901884622323754,
        0.009276957969709997,
    )
    if not np.allclose(q32, expected_q32, rtol=3e-12, atol=3e-14):
        raise AssertionError(("inherited double-cubic formulas", q32))
    discrete = (
        result.closed_entries,
        result.closed_orbits,
        result.extreme_entries,
        result.balanced_entries,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
        result.remaining_quintic_entries,
        result.remaining_extreme_entries,
        result.remaining_balanced_entries,
        len(remaining_quintic_entries()),
    )
    if discrete != (32, 8, 24, 8, 348, 380, 508, 48, 32, 16, 48):
        raise AssertionError(("adjacent double-cubic discrete", discrete))
    observed = (
        result.forward_double_cubic_coefficient,
        result.complement_double_cubic_coefficient,
        result.double_cubic_coefficient,
        result.quintic_fixed_four_energy,
        result.quintic_endpoint_factor,
        result.coefficient,
        result.middle_whole_extreme_coefficient,
        result.quintic_fixed_three_energy,
        result.middle_whole_balanced_coefficient,
        result.generic_extreme_mask_factor,
        result.generic_balanced_mask_factor,
        result.minimum_coefficient,
        result.maximum_coefficient,
        result.previous_routing.total,
        result.inserted_routing.total,
        result.inserted_routing.margin_to_one_third,
        result.routing_change,
        result.remaining_quintic_local_proxy.total,
        result.proxy_reserve_after_declared_allowance,
    )
    expected = (
        7.35015456248367e-07,
        4.704098919989549e-05,
        4.704098919989549e-05,
        0.9990234375,
        0.9995115994824673,
        4.701801435642501e-05,
        0.0002680111923367997,
        1.4538457961309523,
        0.0031270634402608675,
        2.9978026001826406,
        5.1815405503520555,
        4.701801435642501e-05,
        0.013438055203303103,
        0.32734300501757363,
        0.32605080644648143,
        0.007282526886851881,
        -0.0012921985710921957,
        0.3262031888677698,
        0.006130144465563537,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-12):
        raise AssertionError(("adjacent double-cubic numeric", observed))
    committed = (
        ROOT
        / "artifacts"
        / "q64_adjacent_double_cubic_quintic_endpoint_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale adjacent double-cubic endpoint artifact")
    if set(coefficient_map()) != set(entries):
        raise AssertionError("adjacent double-cubic coefficient map")
    target_values = tuple(target_coefficient(target) for target in TARGETS)
    expected_target_values = (
        4.701801435642501e-05,
        0.0002680111923367997,
        0.0031270634402608675,
        0.0013708671889487963,
        0.00014101959973861021,
        0.013438055203303103,
        0.0008038372437911877,
        0.0013708671889487963,
    )
    if not np.allclose(
        target_values, expected_target_values, rtol=3e-9, atol=3e-12
    ):
        raise AssertionError(("adjacent double-cubic target map", target_values))
    print(
        "q64 adjacent double-cubic quintic-endpoint insertion passed: "
        f"entries={result.closed_entries},"
        f"coefficients={result.coefficient:.12g}/"
        f"{result.middle_whole_extreme_coefficient:.12g}/"
        f"{result.middle_whole_balanced_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.inserted_routing.total:.12g},"
        f"margin={result.inserted_routing.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
