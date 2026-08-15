#!/usr/bin/env python3
"""Live quintic gate after all q64 universal coefficient-one insertions."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from math import sqrt
from pathlib import Path

from scipy.optimize import brentq

from q64_paper_target_gate import (
    DIMENSION,
    MODES,
    ORDER,
    RESERVE_TARGET,
    THRESHOLD,
    OptimizedLedger,
    ProfileSplit,
    optimize,
)
import occupation_compatible_sector_optimization as occupation
from q64_remaining_class_gates import contraction_class, remaining_entries
from q64_universal_multicubic_insertion import inserted_coefficients


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class PostUniversalQuinticGate:
    order: int
    dimension: int
    sign_modes: int
    proved_entries: int
    remaining_open_entries: int
    quintic_entries: int
    extreme_split_entries: int
    balanced_split_entries: int
    current_routing: OptimizedLedger
    common_reserve_coefficient: float
    extreme_local_slice_coefficient: float
    balanced_local_slice_coefficient: float
    balanced_reserve_coefficient_after_extreme_local_slice: float
    local_slice_proxy: OptimizedLedger
    local_slice_proxy_overshoot: float


def quintic_entries() -> tuple[ProfileSplit, ...]:
    return tuple(
        entry
        for entry in remaining_entries()
        if contraction_class(entry) == "one_split_cubic_one_split_higher"
        and 5 in entry[0]
    )


def quintic_split_depth(entry: ProfileSplit) -> int:
    profile, split = entry
    index = profile.index(5)
    return min(split[index], 5 - split[index])


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

    return float(brentq(residual, 0.0, 0.3, xtol=2e-10))


def diagnostic() -> PostUniversalQuinticGate:
    entries = quintic_entries()
    extreme = tuple(entry for entry in entries if quintic_split_depth(entry) == 1)
    balanced = tuple(entry for entry in entries if quintic_split_depth(entry) == 2)
    base = inserted_coefficients()
    current = optimize(mapped_coefficients=base)

    cubic_slice = occupation.endpoint_singleton_slice_energies(ORDER)[2]
    quintic_slices = occupation.endpoint_quintic_singleton_slice_energies(ORDER)
    extreme_local = sqrt(cubic_slice * quintic_slices[4])
    balanced_local = sqrt(cubic_slice * quintic_slices[3])

    common_reserve = coefficient_gate(
        base, entries, THRESHOLD - RESERVE_TARGET
    )
    after_extreme = dict(base)
    for entry in extreme:
        after_extreme[entry] = extreme_local
    balanced_reserve = coefficient_gate(
        after_extreme, balanced, THRESHOLD - RESERVE_TARGET
    )

    local_proxy_coefficients = dict(after_extreme)
    for entry in balanced:
        local_proxy_coefficients[entry] = balanced_local
    local_proxy = optimize(mapped_coefficients=local_proxy_coefficients)
    return PostUniversalQuinticGate(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        proved_entries=220,
        remaining_open_entries=668,
        quintic_entries=len(entries),
        extreme_split_entries=len(extreme),
        balanced_split_entries=len(balanced),
        current_routing=current,
        common_reserve_coefficient=common_reserve,
        extreme_local_slice_coefficient=extreme_local,
        balanced_local_slice_coefficient=balanced_local,
        balanced_reserve_coefficient_after_extreme_local_slice=balanced_reserve,
        local_slice_proxy=local_proxy,
        local_slice_proxy_overshoot=-local_proxy.margin_to_one_third,
    )


def artifact_text(result: PostUniversalQuinticGate) -> str:
    payload = {
        "schema": "round4_q64_post_universal_quintic_gate_v1",
        "result": asdict(result),
        "evidence_label": (
            "floating Perron and attenuation optimization after 220 proved "
            "q64 entries; local slice values are unproved all-placement "
            "proxies; one batch only"
        ),
        "decision": (
            "the extreme 1|4 local scale fits, but the balanced 2|3 local "
            "scale misses its live reserve gate; prove additional shared "
            "cancellation or recover margin from a sharper universal insertion"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic post-universal quintic gate",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_post_universal_quintic_gate.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 post-universal quintic gate: "
        f"entries={result.quintic_entries},"
        f"common_reserve={result.common_reserve_coefficient:.12g},"
        f"extreme_local={result.extreme_local_slice_coefficient:.12g},"
        f"balanced_local={result.balanced_local_slice_coefficient:.12g},"
        "balanced_reserve="
        f"{result.balanced_reserve_coefficient_after_extreme_local_slice:.12g},"
        f"proxy_total={result.local_slice_proxy.total:.12g},"
        f"proxy_overshoot={result.local_slice_proxy_overshoot:.12g},"
        "status=routing_gate_not_theorem"
    )


if __name__ == "__main__":
    main()
