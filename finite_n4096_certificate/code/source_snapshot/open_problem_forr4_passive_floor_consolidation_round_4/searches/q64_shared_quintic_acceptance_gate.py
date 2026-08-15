#!/usr/bin/env python3
"""Acceptance gate for the remaining q64 quintics and adaptive interface."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from math import sqrt
from pathlib import Path

from scipy.optimize import brentq

from q64_adjacent_double_cubic_quintic_endpoint_insertion import (
    inserted_coefficients,
    remaining_quintic_entries,
    remaining_quintic_local_proxy_coefficients,
)
from q64_paper_target_gate import (
    ORDER,
    RESERVE_TARGET,
    THRESHOLD,
    OptimizedLedger,
    ProfileSplit,
    optimize,
)
from q64_post_universal_quintic_gate import quintic_split_depth
import occupation_compatible_sector_optimization as occupation


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class SharedQuinticAcceptanceGate:
    order: int
    dimension: int
    proved_entries: int
    remaining_open_entries: int
    remaining_quintic_entries: int
    extreme_quintic_entries: int
    balanced_quintic_entries: int
    remaining_nonquintic_entries: int
    higher_split_only_in_cubic_profile_entries: int
    noncubic_profile_entries: int
    two_split_cubics_one_split_higher_entries: int
    one_split_cubic_no_split_higher_entries: int
    current_routing: OptimizedLedger
    common_reserve_coefficient: float
    extreme_sufficient_coefficient: float
    balanced_sufficient_coefficient: float
    balanced_reserve_gate_after_extreme: float
    sufficient_two_tier_proxy: OptimizedLedger
    raw_adaptive_overhead_cap: float
    adaptive_overhead_cap_retaining_allowance: float
    raw_adaptive_multiplier_cap: float
    adaptive_multiplier_cap_retaining_allowance: float


def coefficient_gate(
    base: dict[ProfileSplit, float],
    entries: tuple[ProfileSplit, ...],
    target_total: float,
) -> float:
    cache: dict[float, float] = {}

    def residual(value: float) -> float:
        if value not in cache:
            trial = dict(base)
            for entry in entries:
                trial[entry] = value
            cache[value] = optimize(mapped_coefficients=trial).total
        return cache[value] - target_total

    upper = 0.3
    while residual(upper) < 0 and upper < 4.8:
        upper *= 2
    if residual(upper) < 0:
        raise ValueError(("coefficient gate exceeds search bracket", upper))
    return float(brentq(residual, 0.0, upper, xtol=2e-10))


def diagnostic() -> SharedQuinticAcceptanceGate:
    entries = remaining_quintic_entries()
    extreme_entries = tuple(
        entry for entry in entries if quintic_split_depth(entry) == 1
    )
    balanced_entries = tuple(
        entry for entry in entries if quintic_split_depth(entry) == 2
    )
    base = inserted_coefficients()
    current = optimize(mapped_coefficients=base)
    cubic_slice = occupation.endpoint_singleton_slice_energies(ORDER)[2]
    quintic_slices = occupation.endpoint_quintic_singleton_slice_energies(
        ORDER
    )
    extreme = sqrt(cubic_slice * quintic_slices[4])
    balanced = sqrt(cubic_slice * quintic_slices[3])
    common_gate = coefficient_gate(
        base, entries, THRESHOLD - RESERVE_TARGET
    )
    after_extreme = dict(base)
    for entry in extreme_entries:
        after_extreme[entry] = extreme
    balanced_gate = coefficient_gate(
        after_extreme,
        balanced_entries,
        THRESHOLD - RESERVE_TARGET,
    )
    local_proxy = optimize(
        mapped_coefficients=remaining_quintic_local_proxy_coefficients()
    )
    raw_overhead = local_proxy.margin_to_one_third
    retained_overhead = raw_overhead - RESERVE_TARGET
    return SharedQuinticAcceptanceGate(
        order=ORDER,
        dimension=ORDER * ORDER,
        proved_entries=380,
        remaining_open_entries=508,
        remaining_quintic_entries=len(entries),
        extreme_quintic_entries=len(extreme_entries),
        balanced_quintic_entries=len(balanced_entries),
        remaining_nonquintic_entries=460,
        higher_split_only_in_cubic_profile_entries=176,
        noncubic_profile_entries=140,
        two_split_cubics_one_split_higher_entries=96,
        one_split_cubic_no_split_higher_entries=48,
        current_routing=current,
        common_reserve_coefficient=common_gate,
        extreme_sufficient_coefficient=extreme,
        balanced_sufficient_coefficient=balanced,
        balanced_reserve_gate_after_extreme=balanced_gate,
        sufficient_two_tier_proxy=local_proxy,
        raw_adaptive_overhead_cap=raw_overhead,
        adaptive_overhead_cap_retaining_allowance=retained_overhead,
        raw_adaptive_multiplier_cap=THRESHOLD / local_proxy.total,
        adaptive_multiplier_cap_retaining_allowance=(
            (THRESHOLD - RESERVE_TARGET) / local_proxy.total
        ),
    )


def artifact_text(result: SharedQuinticAcceptanceGate) -> str:
    payload = {
        "schema": "round4_q64_shared_quintic_acceptance_gate_v1",
        "result": asdict(result),
        "evidence_label": (
            "floating q64 Perron routing gates for the 48 remaining "
            "split-cubic/split-quintic entries; adaptive caps are algebraic "
            "requirements, not an adaptive theorem"
        ),
        "acceptance": (
            "prove arbitrary-law coverage of all 48 entries and an "
            "outward-rounded routing total at most 1/3-1e-3 with explicit "
            "gates for the other 460 entries; derive an "
            "outcome-width/depth-uniform adaptive recurrence whose evaluated "
            "overhead fits inside the certified remaining margin"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic shared-quintic acceptance gate",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_shared_quintic_acceptance_gate.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 shared quintic acceptance gate: "
        f"entries={result.remaining_quintic_entries},"
        f"extreme={result.extreme_quintic_entries},"
        f"balanced={result.balanced_quintic_entries},"
        f"common_reserve={result.common_reserve_coefficient:.12g},"
        "balanced_reserve_after_extreme="
        f"{result.balanced_reserve_gate_after_extreme:.12g},"
        f"two_tier_total={result.sufficient_two_tier_proxy.total:.12g},"
        "adaptive_additive_cap="
        f"{result.adaptive_overhead_cap_retaining_allowance:.12g},"
        "adaptive_multiplier_cap="
        f"{result.adaptive_multiplier_cap_retaining_allowance:.12g},"
        "status=routing_and_adaptive_requirement_not_theorem"
    )


if __name__ == "__main__":
    main()
