#!/usr/bin/env python3
"""Regression for the q64 shared-quintic row/chain insertion."""

from __future__ import annotations

from itertools import combinations
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

from adjacent_cubic_quintic_orbit_witness import (  # noqa: E402
    parity_record_size,
    record_one_link_moment,
)
from q64_adjacent_double_cubic_quintic_endpoint_insertion import (  # noqa: E402
    remaining_quintic_entries as pre_insertion_entries,
)
from q64_degree_ten_completion_row_insertion import orbit  # noqa: E402
from q64_paper_target_gate import (  # noqa: E402
    RESERVE_TARGET,
    THRESHOLD,
    optimize,
)
from q64_remaining_class_gates import partition_remaining  # noqa: E402
from q64_shared_quintic_row_chain_insertion import (  # noqa: E402
    CUBIC_ENDPOINT_FIXED_FOUR,
    FAVORABLE_ADJACENT_SPLIT,
    FAVORABLE_ADJACENT_WHOLE,
    FIXED_ONE_ADJACENT_WHOLE,
    NONFAVORABLE_ADJACENT_SPLIT,
    RESIDUAL_CLASS_LABELS,
    artifact_text,
    coefficient_map,
    diagnostic,
    endpoint_quintic_fixed_one_record_energies,
    odd_record_incidence,
    remaining_quintic_entries,
    shared_quintic_entries,
    inserted_coefficients,
)


def direct_q4_fixed_one_record_energies() -> tuple[float, float, float]:
    """Enumerate the fixed-one endpoint slice by its other record."""

    order = 4
    dimension = order * order
    energies = {1: 0.0, 3: 0.0, 5: 0.0}
    for remainder in combinations(range(1, dimension), 4):
        quintic = (0,) + remainder
        moment = record_one_link_moment(order, (0,), quintic)
        other_record = parity_record_size(order, quintic, axis=1)
        energies[other_record] += moment * moment
    return tuple(energies[record] for record in (1, 3, 5))


def main() -> None:
    result = diagnostic()
    entries = shared_quintic_entries()
    if set(entries) != set(pre_insertion_entries()):
        raise AssertionError("shared-quintic family does not close exact remainder")
    if remaining_quintic_entries():
        raise AssertionError(("unexpected remaining quintics", remaining_quintic_entries()))
    family_counts = tuple(
        len({entry for target in family for entry in orbit(target)})
        for family in (
            FAVORABLE_ADJACENT_SPLIT,
            FAVORABLE_ADJACENT_WHOLE,
            NONFAVORABLE_ADJACENT_SPLIT,
            FIXED_ONE_ADJACENT_WHOLE,
            CUBIC_ENDPOINT_FIXED_FOUR,
        )
    )
    if family_counts != (16, 12, 8, 8, 4):
        raise AssertionError(("shared-quintic family counts", family_counts))
    q4_formula = endpoint_quintic_fixed_one_record_energies(4)
    q4_direct = direct_q4_fixed_one_record_energies()
    if not np.allclose(q4_formula, q4_direct, rtol=3e-13, atol=3e-13):
        raise AssertionError(("q4 record-energy identity", q4_formula, q4_direct))
    if not np.isclose(sum(q4_formula), 85 / 16, rtol=0, atol=3e-13):
        raise AssertionError(("q4 fixed-one total", q4_formula))
    if (
        odd_record_incidence(64, 5, 1, 4),
        odd_record_incidence(64, 5, 3, 4),
    ) != (4092, 3968):
        raise AssertionError("fixed-four record incidences")
    discrete = (
        result.closed_entries,
        result.closed_orbits,
        result.extreme_entries,
        result.balanced_entries,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
        result.remaining_quintic_entries,
    )
    if discrete != (48, 12, 32, 16, 380, 428, 460, 0):
        raise AssertionError(("shared-quintic discrete", discrete))
    observed = (
        result.fixed_one_record_one_energy,
        result.fixed_one_record_three_energy,
        result.fixed_one_record_five_energy,
        result.fixed_one_total_energy,
        result.record_one_middle_maximum,
        result.record_three_middle_maximum,
        result.favorable_adjacent_split_extreme_coefficient,
        result.favorable_adjacent_split_balanced_coefficient,
        result.favorable_adjacent_whole_extreme_coefficient,
        result.favorable_adjacent_whole_balanced_coefficient,
        result.nonfavorable_adjacent_split_extreme_coefficient,
        result.nonfavorable_adjacent_split_balanced_coefficient,
        result.fixed_one_adjacent_whole_coefficient,
        result.cubic_endpoint_fixed_four_coefficient,
        result.minimum_coefficient,
        result.maximum_coefficient,
        result.previous_routing.total,
        result.inserted_routing.total,
        result.inserted_routing.margin_to_one_third,
        result.routing_improvement,
        result.adaptive_additive_cap_retaining_allowance,
        result.adaptive_multiplier_cap_retaining_allowance,
    )
    expected = (
        1.2536621093750002,
        55.77392578125,
        968.22265625,
        1025.250244140625,
        0.0002640168970814132,
        1.180403414670998e-06,
        0.018927574738987724,
        0.020373745136755918,
        0.0002638879510922413,
        0.00031833976776181127,
        0.0183664518516331,
        0.020210337661283738,
        0.008453696720453088,
        5.350510859088515e-06,
        5.350510859088515e-06,
        0.020373745136755918,
        0.32605080644648143,
        0.3238115631713356,
        0.009521770161997734,
        0.0022392432751458524,
        0.008521770161997733,
        1.0263170656369942,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-12):
        raise AssertionError(("shared-quintic numeric", observed))
    if set(coefficient_map()) != set(entries):
        raise AssertionError("shared-quintic coefficient map")
    partition = partition_remaining()
    base = inserted_coefficients()
    for label, count, frozen, gate in zip(
        result.residual_class_labels,
        result.residual_class_counts,
        result.residual_class_frozen_targets,
        result.residual_class_reserve_gates,
        strict=True,
    ):
        if label not in RESIDUAL_CLASS_LABELS:
            raise AssertionError(("residual class label", label))
        class_entries = partition[label]
        if len(class_entries) != count:
            raise AssertionError(("residual class count", label, count))
        if not all(base[entry] == frozen for entry in class_entries):
            raise AssertionError(("residual frozen target", label, frozen))
        trial = dict(base)
        for entry in class_entries:
            trial[entry] = gate
        gate_total = optimize(mapped_coefficients=trial).total
        if not np.isclose(
            gate_total,
            THRESHOLD - RESERVE_TARGET,
            rtol=0,
            atol=4e-10,
        ):
            raise AssertionError(("residual reserve gate", label, gate_total))
    committed = (
        ROOT / "artifacts" / "q64_shared_quintic_row_chain_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale shared-quintic row/chain artifact")
    print(
        "q64 shared quintic row/chain insertion passed: "
        f"entries={result.closed_entries},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.inserted_routing.total:.12g},"
        f"margin={result.inserted_routing.margin_to_one_third:.12g},"
        "adaptive_additive_cap="
        f"{result.adaptive_additive_cap_retaining_allowance:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
