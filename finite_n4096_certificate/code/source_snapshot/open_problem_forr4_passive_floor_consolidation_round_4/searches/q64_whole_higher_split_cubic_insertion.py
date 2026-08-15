#!/usr/bin/env python3
"""Close 48 q64 entries with one split cubic and whole higher blocks."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from json import dumps
from math import sqrt
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent
    / "open_problem_forr4_passive_floor_consolidation_round_3"
    / "searches"
)
sys.path.insert(0, str(ROUND3_SEARCHES))

from middle_cubic_quintic_pair_contraction import (  # noqa: E402
    middle_link_maxima,
)
import occupation_compatible_sector_optimization as occupation  # noqa: E402
from q64_dual_endpoint_schur_insertion import (  # noqa: E402
    has_favorable_cubic_singleton,
)
from q64_noncubic_recovered_universal_insertion import (  # noqa: E402
    inserted_coefficients as noncubic_inserted_coefficients,
)
from q64_paper_target_gate import (  # noqa: E402
    DIMENSION,
    MODES,
    ORDER,
    RESERVE_TARGET,
    THRESHOLD,
    OptimizedLedger,
    ProfileSplit,
    optimize,
)
from q64_remaining_class_gates import partition_remaining  # noqa: E402
from q64_shared_quintic_row_chain_insertion import (  # noqa: E402
    endpoint_quintic_fixed_one_record_energies,
)


CLASS_LABEL = "one_split_cubic_no_split_higher"
REMAINING_CLASS_LABEL = "higher_split_only_in_cubic_profile"
REMAINING_CLASS_RESERVE_GATE = 0.1425819092113566


@dataclass(frozen=True)
class WholeHigherSplitCubicInsertion:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    favorable_endpoint_entries: int
    internal_endpoint_entries: int
    complete_wedge_entries: int
    favorable_endpoint_coefficient: float
    internal_endpoint_coefficient: float
    endpoint_record_one_full_energy: float
    endpoint_record_three_full_energy: float
    endpoint_record_five_full_energy: float
    record_one_middle_maximum: float
    record_three_middle_maximum: float
    complete_wedge_coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    inserted_routing: OptimizedLedger
    routing_improvement: float
    reserve_after_declared_allowance: float
    adaptive_multiplier_cap_retaining_allowance: float
    remaining_class_entries: int
    remaining_class_frozen_target: float
    remaining_class_reserve_gate: float


def target_entries() -> tuple[ProfileSplit, ...]:
    return partition_remaining()[CLASS_LABEL]


def has_internal_singleton_cubic(entry: ProfileSplit) -> bool:
    """Whether a whole cubic shares the singleton side on their link."""

    profile, split = entry
    singleton = profile.index(1)
    side = split[singleton]
    return any(
        profile[neighbor] == 3
        and split[neighbor] in (0, 3)
        and split[neighbor] == 3 * side
        for neighbor in (singleton - 1, singleton + 1)
        if 0 <= neighbor < 4
    )


def entry_type(entry: ProfileSplit) -> str:
    if has_favorable_cubic_singleton(entry):
        return "favorable_endpoint"
    if has_internal_singleton_cubic(entry):
        return "internal_endpoint"
    return "complete_wedge"


def full_endpoint_record_energies(
    order: int = ORDER,
) -> tuple[float, float, float]:
    """Convert fixed-one energies to full endpoint row energies."""

    fixed_one = endpoint_quintic_fixed_one_record_energies(order)
    return tuple(value * order**2 / 5 for value in fixed_one)


def favorable_endpoint_coefficient(order: int = ORDER) -> float:
    energy = occupation.endpoint_singleton_slice_energies(order)[2]
    return sqrt(energy)


def complete_wedge_coefficient(order: int = ORDER) -> float:
    """Weighted complete-row coefficient for ``M_15 M_53``."""

    record_energies = full_endpoint_record_energies(order)
    record_one, record_three, _ = middle_link_maxima(order)
    row_energy = (
        record_energies[0] * record_one**2
        + record_energies[1] * record_three**2
    )
    # The record-five endpoint sector cannot couple to a cubic. The universal
    # cross-Gram theorem remains a safe fallback at small orders.
    return min(1.0, sqrt(row_energy))


def coefficient_map() -> dict[ProfileSplit, float]:
    values = {
        "favorable_endpoint": favorable_endpoint_coefficient(),
        "internal_endpoint": 1 / ORDER,
        "complete_wedge": complete_wedge_coefficient(),
    }
    return {entry: values[entry_type(entry)] for entry in target_entries()}


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = noncubic_inserted_coefficients()
    result.update(coefficient_map())
    return result


def diagnostic() -> WholeHigherSplitCubicInsertion:
    entries = target_entries()
    types = Counter(entry_type(entry) for entry in entries)
    full_energies = full_endpoint_record_energies()
    record_one, record_three, _ = middle_link_maxima(ORDER)
    previous = optimize(mapped_coefficients=noncubic_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    remaining = partition_remaining()[REMAINING_CLASS_LABEL]
    inserted_map = inserted_coefficients()
    previous_proved = 664
    return WholeHigherSplitCubicInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        favorable_endpoint_entries=types["favorable_endpoint"],
        internal_endpoint_entries=types["internal_endpoint"],
        complete_wedge_entries=types["complete_wedge"],
        favorable_endpoint_coefficient=favorable_endpoint_coefficient(),
        internal_endpoint_coefficient=1 / ORDER,
        endpoint_record_one_full_energy=full_energies[0],
        endpoint_record_three_full_energy=full_energies[1],
        endpoint_record_five_full_energy=full_energies[2],
        record_one_middle_maximum=record_one,
        record_three_middle_maximum=record_three,
        complete_wedge_coefficient=complete_wedge_coefficient(),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        inserted_routing=inserted,
        routing_improvement=previous.total - inserted.total,
        reserve_after_declared_allowance=(
            inserted.margin_to_one_third - RESERVE_TARGET
        ),
        adaptive_multiplier_cap_retaining_allowance=(
            (THRESHOLD - RESERVE_TARGET) / inserted.total
        ),
        remaining_class_entries=len(remaining),
        remaining_class_frozen_target=inserted_map[remaining[0]],
        remaining_class_reserve_gate=REMAINING_CLASS_RESERVE_GATE,
    )


def artifact_text(result: WholeHigherSplitCubicInsertion) -> str:
    payload = {
        "schema": "round4_q64_whole_higher_split_cubic_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal q64 theorem for 48 entries with one split "
            "cubic and whole higher blocks: exact cubic endpoint slice, "
            "internal singleton-cubic scalar, or complete M15-M53 weighted "
            "row followed by unit cross-Gram dressing; floating Perron "
            "insertion; one batch only"
        ),
        "adaptive_requirement": (
            "the local 48-entry theorem survives; its cumulative count, "
            "routing reserve, and multiplier are withdrawn because the "
            "upstream universal insertion is quarantined"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic whole-higher split-cubic artifact",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = (
            ROOT
            / "artifacts"
            / "q64_whole_higher_split_cubic_insertion.json"
        )
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 whole-higher split-cubic insertion: "
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
