#!/usr/bin/env python3
"""Compute the lead common-coefficient gate after the q64 insertions.

The q64 block-coherent and inherited chain-aware theorems close 110 of the
888 balanced high-sector entries.  This module partitions the other 778 by
the number of internally split cubic and higher-degree blocks.  For each
class it records the exact inventory and existing routing target.  For the
dominant 280-entry class it varies one common coefficient, leaves every other
open entry at the existing target, and reoptimizes the occupation ledger.

The resulting values are proof-allocation gates, not theorem coefficients.
They answer how loose the next shared contraction may be while retaining
either zero or 1e-3 reserve.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from pathlib import Path

from scipy.optimize import brentq

from q64_block_coherent_contraction import block_coherent_entries
from q64_chain_aware_insertion import (
    chain_aware_entries,
    inserted_coefficients,
)
from q64_paper_target_gate import (
    DIMENSION,
    MODES,
    ORDER,
    RESERVE_TARGET,
    THRESHOLD,
    OptimizedLedger,
    ProfileSplit,
    balanced_open_entries,
    optimize,
)


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class RemainingClassRow:
    name: str
    split_cubic_blocks_min: int
    split_cubic_blocks_max: int
    split_higher_blocks_min: int
    split_higher_blocks_max: int
    entries: int
    routing_target: float


@dataclass(frozen=True)
class RemainingClassGates:
    order: int
    dimension: int
    sign_modes: int
    balanced_open_entries: int
    proved_entries: int
    remaining_entries: int
    baseline: OptimizedLedger
    classes: tuple[RemainingClassRow, ...]
    lead_class: str
    lead_entries: int
    lead_threshold_coefficient: float
    lead_reserve_coefficient: float
    lead_reserve_multiplier_over_target: float


CLASS_ORDER = (
    "one_split_cubic_one_split_higher",
    "higher_split_only_in_cubic_profile",
    "noncubic_profile",
    "two_split_cubics_one_split_higher",
    "one_split_cubic_no_split_higher",
    "two_split_cubics_no_split_higher",
    "three_split_cubics",
    "four_split_cubics",
)

COUNT_CLASS_NAMES = {
    (1, 1): "one_split_cubic_one_split_higher",
    (0, 1): "higher_split_only_in_cubic_profile",
    (2, 1): "two_split_cubics_one_split_higher",
    (1, 0): "one_split_cubic_no_split_higher",
    (2, 0): "two_split_cubics_no_split_higher",
    (3, 0): "three_split_cubics",
    (4, 0): "four_split_cubics",
}


def internal_split_counts(entry: ProfileSplit) -> tuple[int, int]:
    profile, split = entry
    cubic = sum(
        degree == 3 and selected not in (0, degree)
        for degree, selected in zip(profile, split, strict=True)
    )
    higher = sum(
        degree > 3 and selected not in (0, degree)
        for degree, selected in zip(profile, split, strict=True)
    )
    return cubic, higher


def remaining_entries() -> tuple[ProfileSplit, ...]:
    proved = set(block_coherent_entries()).union(chain_aware_entries())
    return tuple(entry for entry in balanced_open_entries() if entry not in proved)


def contraction_class(entry: ProfileSplit) -> str:
    profile, _ = entry
    if 3 not in profile:
        return "noncubic_profile"
    counts = internal_split_counts(entry)
    if counts not in COUNT_CLASS_NAMES:
        raise AssertionError(("unexpected cubic-profile class", entry, counts))
    return COUNT_CLASS_NAMES[counts]


def partition_remaining() -> dict[str, tuple[ProfileSplit, ...]]:
    partition: dict[str, list[ProfileSplit]] = {}
    for entry in remaining_entries():
        label = contraction_class(entry)
        partition.setdefault(label, []).append(entry)
    if set(partition) != set(CLASS_ORDER):
        raise AssertionError(("unexpected remaining classes", sorted(partition)))
    return {label: tuple(value) for label, value in partition.items()}


def class_gate(
    entries: tuple[ProfileSplit, ...], target_total: float
) -> float:
    base = inserted_coefficients()
    cache: dict[float, float] = {}

    def residual(value: float) -> float:
        if value not in cache:
            trial = dict(base)
            for entry in entries:
                trial[entry] = value
            cache[value] = optimize(mapped_coefficients=trial).total
        return cache[value] - target_total

    if residual(0.0) >= 0:
        raise AssertionError(("class cannot pass even at zero", len(entries)))
    upper = 0.25
    while residual(upper) < 0:
        upper *= 2
    return float(brentq(residual, 0.0, upper, xtol=2e-10))


def diagnostic() -> RemainingClassGates:
    partition = partition_remaining()
    base = inserted_coefficients()
    baseline = optimize(mapped_coefficients=base)
    rows = []
    for label in CLASS_ORDER:
        entries = partition[label]
        counts = tuple(internal_split_counts(entry) for entry in entries)
        routing_target = base[entries[0]]
        if not all(base[entry] == routing_target for entry in entries):
            raise AssertionError(("nonuniform routing target", label))
        rows.append(
            RemainingClassRow(
                name=label,
                split_cubic_blocks_min=min(value[0] for value in counts),
                split_cubic_blocks_max=max(value[0] for value in counts),
                split_higher_blocks_min=min(value[1] for value in counts),
                split_higher_blocks_max=max(value[1] for value in counts),
                entries=len(entries),
                routing_target=routing_target,
            )
        )
    lead_label = CLASS_ORDER[0]
    lead_entries = partition[lead_label]
    lead_target = base[lead_entries[0]]
    threshold = class_gate(lead_entries, THRESHOLD)
    reserve = class_gate(lead_entries, THRESHOLD - RESERVE_TARGET)
    return RemainingClassGates(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        balanced_open_entries=len(balanced_open_entries()),
        proved_entries=len(block_coherent_entries()) + len(chain_aware_entries()),
        remaining_entries=len(remaining_entries()),
        baseline=baseline,
        classes=tuple(rows),
        lead_class=lead_label,
        lead_entries=len(lead_entries),
        lead_threshold_coefficient=threshold,
        lead_reserve_coefficient=reserve,
        lead_reserve_multiplier_over_target=reserve / lead_target,
    )


def artifact_text(result: RemainingClassGates) -> str:
    payload = {
        "schema": "round4_q64_remaining_class_gates_v1",
        "result": asdict(result),
        "evidence_label": (
            "floating Perron and attenuation optimization after 110 proved "
            "q64 entries; displayed class coefficients are routing gates, "
            "not arbitrary-law theorems; one batch only"
        ),
        "next_target": (
            "prove or falsify one shared contraction for the 280-entry "
            "one-split-cubic/one-split-higher class at its reserve gate"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write deterministic remaining-class gates under artifacts/",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_remaining_class_gates.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 remaining class gates: "
        f"proved={result.proved_entries},"
        f"remaining={result.remaining_entries},"
        f"baseline={result.baseline.total:.12g},"
        f"baseline_margin={result.baseline.margin_to_one_third:.12g},"
        f"lead_entries={result.lead_entries},"
        f"lead_threshold={result.lead_threshold_coefficient:.12g},"
        f"lead_reserve={result.lead_reserve_coefficient:.12g},"
        f"lead_multiplier={result.lead_reserve_multiplier_over_target:.12g},"
        "status=routing_gates_not_theorems"
    )


if __name__ == "__main__":
    main()
