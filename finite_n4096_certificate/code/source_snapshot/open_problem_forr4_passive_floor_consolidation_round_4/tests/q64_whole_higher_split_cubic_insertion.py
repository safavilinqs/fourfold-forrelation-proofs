#!/usr/bin/env python3
"""Regression for the q64 whole-higher split-cubic insertion."""

from __future__ import annotations

from itertools import combinations, permutations, product
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent
    / "open_problem_forr4_passive_floor_consolidation_round_3"
    / "searches"
)
sys.path.insert(0, str(ROOT / "searches"))
sys.path.insert(0, str(ROUND3_SEARCHES))

from q64_paper_target_gate import (  # noqa: E402
    RESERVE_TARGET,
    THRESHOLD,
    optimize,
)
from q64_remaining_class_gates import partition_remaining  # noqa: E402
from q64_whole_higher_split_cubic_insertion import (  # noqa: E402
    REMAINING_CLASS_LABEL,
    artifact_text,
    coefficient_map,
    complete_wedge_coefficient,
    diagnostic,
    entry_type,
    full_endpoint_record_energies,
    inserted_coefficients,
    target_entries,
)


def direct_q4_wedge_energy() -> float:
    """Exact max sum_Q |M_15(S,Q) M_53(Q,C)|^2 at q=4."""

    order = 4
    dimension = order * order
    hadamard = np.asarray([[1]], dtype=np.int8)
    while len(hadamard) < order:
        hadamard = np.block(
            [[hadamard, hadamard], [hadamard, -hadamard]]
        )
    left = []
    right = []
    for permutation in permutations(range(order)):
        for signs in product((-1, 1), repeat=order):
            signed_permutation = np.zeros(
                (order, order), dtype=np.int8
            )
            for column, row in enumerate(permutation):
                signed_permutation[row, column] = signs[column]
            left.append((hadamard @ signed_permutation).reshape(-1))
            right.append((signed_permutation @ hadamard).reshape(-1))
    left_array = np.asarray(left, dtype=np.int8)
    right_array = np.asarray(right, dtype=np.int8)
    cubics = np.asarray(
        tuple(combinations(range(dimension), 3)), dtype=np.int16
    )
    quintics = np.asarray(
        tuple(combinations(range(dimension), 5)), dtype=np.int16
    )
    left_five = np.prod(
        left_array[:, quintics], axis=2, dtype=np.int8
    )
    right_three = np.prod(
        right_array[:, cubics], axis=2, dtype=np.int8
    )
    right_five = np.prod(
        right_array[:, quintics], axis=2, dtype=np.int8
    )
    normalization = len(left_array)
    moment_15 = left_array.T.astype(float) @ right_five / normalization
    moment_53 = left_five.T.astype(float) @ right_three / normalization
    energy = (moment_15 * moment_15) @ (moment_53 * moment_53)
    return float(np.max(energy))


def main() -> None:
    result = diagnostic()
    entries = target_entries()
    if set(entries) != set(
        partition_remaining()["one_split_cubic_no_split_higher"]
    ):
        raise AssertionError("whole-higher split-cubic inventory")
    types = {label: 0 for label in (
        "favorable_endpoint",
        "internal_endpoint",
        "complete_wedge",
    )}
    for entry in entries:
        types[entry_type(entry)] += 1
    if types != {
        "favorable_endpoint": 24,
        "internal_endpoint": 16,
        "complete_wedge": 8,
    }:
        raise AssertionError(("whole-higher type partition", types))
    coefficient_values = coefficient_map()
    if set(coefficient_values) != set(entries):
        raise AssertionError("whole-higher coefficient coverage")
    if not np.allclose(
        full_endpoint_record_energies(4),
        (7.0, 10.0, 0.0),
        rtol=0,
        atol=3e-14,
    ):
        raise AssertionError("q4 full endpoint record energies")
    direct_wedge = direct_q4_wedge_energy()
    if direct_wedge > complete_wedge_coefficient(4) ** 2 + 3e-14:
        raise AssertionError(("q4 complete wedge", direct_wedge))
    if not np.isclose(
        direct_wedge, 17 / 72, rtol=0, atol=3e-14
    ):
        raise AssertionError(("q4 complete wedge identity", direct_wedge))

    discrete = (
        result.closed_entries,
        result.favorable_endpoint_entries,
        result.internal_endpoint_entries,
        result.complete_wedge_entries,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
        result.remaining_class_entries,
    )
    if discrete != (48, 24, 16, 8, 664, 712, 176, 176):
        raise AssertionError(("whole-higher discrete", discrete))
    observed = (
        result.favorable_endpoint_coefficient,
        result.internal_endpoint_coefficient,
        result.endpoint_record_one_full_energy,
        result.endpoint_record_three_full_energy,
        result.endpoint_record_five_full_energy,
        result.record_one_middle_maximum,
        result.record_three_middle_maximum,
        result.complete_wedge_coefficient,
        result.previous_routing.total,
        result.inserted_routing.total,
        result.inserted_routing.margin_to_one_third,
        result.routing_improvement,
        result.reserve_after_declared_allowance,
        result.adaptive_multiplier_cap_retaining_allowance,
        result.remaining_class_frozen_target,
        result.remaining_class_reserve_gate,
    )
    expected = (
        0.12403521525363623,
        0.015625,
        1027.0000000000002,
        45690.0,
        793168.0,
        0.0002640168970814132,
        1.180403414670998e-06,
        0.008464668753116878,
        0.3289382301229411,
        0.3284774211661729,
        0.004855912167160414,
        0.00046080895676819944,
        0.0038559121671604144,
        1.011738743422519,
        0.12403521525363623,
        0.1425819092113566,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-12):
        raise AssertionError(("whole-higher numeric", observed))

    remaining = partition_remaining()[REMAINING_CLASS_LABEL]
    base = inserted_coefficients()
    if not all(
        base[entry] == result.remaining_class_frozen_target
        for entry in remaining
    ):
        raise AssertionError("remaining class frozen target")
    trial = dict(base)
    for entry in remaining:
        trial[entry] = result.remaining_class_reserve_gate
    gate_total = optimize(mapped_coefficients=trial).total
    if not np.isclose(
        gate_total,
        THRESHOLD - RESERVE_TARGET,
        rtol=0,
        atol=4e-10,
    ):
        raise AssertionError(("remaining class reserve gate", gate_total))

    committed = (
        ROOT
        / "artifacts"
        / "q64_whole_higher_split_cubic_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale whole-higher split-cubic artifact")
    print(
        "q64 whole-higher split-cubic insertion passed: "
        f"entries={result.closed_entries},"
        f"types={result.favorable_endpoint_entries}/"
        f"{result.internal_endpoint_entries}/"
        f"{result.complete_wedge_entries},"
        f"coefficients={result.favorable_endpoint_coefficient:.12g}/"
        f"{result.internal_endpoint_coefficient:.12g}/"
        f"{result.complete_wedge_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.inserted_routing.total:.12g},"
        f"margin={result.inserted_routing.margin_to_one_third:.12g},"
        f"remaining_open={result.remaining_open_entries},"
        "status=local_48_entry_theorem_cumulative_values_withdrawn"
    )


if __name__ == "__main__":
    main()
