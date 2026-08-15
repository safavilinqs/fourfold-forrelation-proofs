#!/usr/bin/env python3
"""Close the 14 lowest-impact multicubic q64 entries by Gram coefficient one."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from pathlib import Path

from q64_paper_target_gate import (
    DIMENSION,
    MODES,
    ORDER,
    RESERVE_TARGET,
    OptimizedLedger,
    ProfileSplit,
    optimize,
)
from q64_remaining_class_gates import partition_remaining
from q64_universal_septimic_insertion import (
    UNIVERSAL_GRAM_COEFFICIENT,
    inserted_coefficients as septimic_inserted_coefficients,
)


ROOT = Path(__file__).resolve().parents[1]
TARGET_CLASSES = ("three_split_cubics", "four_split_cubics")


@dataclass(frozen=True)
class UniversalMulticubicInsertion:
    order: int
    dimension: int
    sign_modes: int
    universal_coefficient: float
    three_split_cubic_entries: int
    four_split_cubic_entries: int
    newly_closed_entries: int
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    multicubic_inserted: OptimizedLedger
    margin_spent: float
    reserve_after_declared_allowance: float


def multicubic_entries() -> tuple[ProfileSplit, ...]:
    partition = partition_remaining()
    return tuple(entry for label in TARGET_CLASSES for entry in partition[label])


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = septimic_inserted_coefficients()
    for entry in multicubic_entries():
        result[entry] = UNIVERSAL_GRAM_COEFFICIENT
    return result


def diagnostic() -> UniversalMulticubicInsertion:
    partition = partition_remaining()
    three = partition["three_split_cubics"]
    four = partition["four_split_cubics"]
    previous = optimize(mapped_coefficients=septimic_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    newly_closed = len(three) + len(four)
    previous_proved = 206
    return UniversalMulticubicInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        universal_coefficient=UNIVERSAL_GRAM_COEFFICIENT,
        three_split_cubic_entries=len(three),
        four_split_cubic_entries=len(four),
        newly_closed_entries=newly_closed,
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + newly_closed,
        remaining_open_entries=888 - previous_proved - newly_closed,
        previous_routing=previous,
        multicubic_inserted=inserted,
        margin_spent=(
            previous.margin_to_one_third - inserted.margin_to_one_third
        ),
        reserve_after_declared_allowance=(
            inserted.margin_to_one_third - RESERVE_TARGET
        ),
    )


def artifact_text(result: UniversalMulticubicInsertion) -> str:
    payload = {
        "schema": "round4_q64_universal_multicubic_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "quarantined unmasked cross-Gram diagnostic for 14 three-/four-"
            "split-cubic entries; it omits cross-cut distinctness masks"
        ),
        "remaining_open": (
            "historical conditional count only; use the masked-universal audit"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic universal multicubic insertion",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_universal_multicubic_insertion.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 universal multicubic insertion: "
        f"entries={result.newly_closed_entries},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.multicubic_inserted.total:.12g},"
        f"margin={result.multicubic_inserted.margin_to_one_third:.12g},"
        "reserve_after_allowance="
        f"{result.reserve_after_declared_allowance:.12g},"
        "status=quarantined_unmasked_coefficient_one_diagnostic"
    )


if __name__ == "__main__":
    main()
