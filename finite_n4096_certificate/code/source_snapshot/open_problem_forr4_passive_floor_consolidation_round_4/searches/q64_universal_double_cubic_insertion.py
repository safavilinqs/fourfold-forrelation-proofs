#!/usr/bin/env python3
"""Close the 24 two-split-cubic q64 entries by Gram coefficient one."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from pathlib import Path

from q64_fixed_singleton_pair_contraction import (
    inserted_coefficients as fixed_singleton_inserted_coefficients,
)
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
from q64_universal_septimic_insertion import UNIVERSAL_GRAM_COEFFICIENT


ROOT = Path(__file__).resolve().parents[1]
TARGET_CLASS = "two_split_cubics_no_split_higher"


@dataclass(frozen=True)
class UniversalDoubleCubicInsertion:
    order: int
    dimension: int
    sign_modes: int
    universal_coefficient: float
    newly_closed_entries: int
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    double_cubic_inserted: OptimizedLedger
    margin_spent: float
    reserve_after_declared_allowance: float


def double_cubic_entries() -> tuple[ProfileSplit, ...]:
    return partition_remaining()[TARGET_CLASS]


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = fixed_singleton_inserted_coefficients()
    for entry in double_cubic_entries():
        result[entry] = UNIVERSAL_GRAM_COEFFICIENT
    return result


def diagnostic() -> UniversalDoubleCubicInsertion:
    entries = double_cubic_entries()
    previous = optimize(mapped_coefficients=fixed_singleton_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    previous_proved = 236
    return UniversalDoubleCubicInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        universal_coefficient=UNIVERSAL_GRAM_COEFFICIENT,
        newly_closed_entries=len(entries),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        double_cubic_inserted=inserted,
        margin_spent=(
            previous.margin_to_one_third - inserted.margin_to_one_third
        ),
        reserve_after_declared_allowance=(
            inserted.margin_to_one_third - RESERVE_TARGET
        ),
    )


def artifact_text(result: UniversalDoubleCubicInsertion) -> str:
    payload = {
        "schema": "round4_q64_universal_double_cubic_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "quarantined unmasked cross-Gram diagnostic for 24 double-cubic "
            "entries; it omits both cross-cut distinctness masks"
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
        help="write the deterministic universal double-cubic insertion",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_universal_double_cubic_insertion.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 universal double-cubic insertion: "
        f"entries={result.newly_closed_entries},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.double_cubic_inserted.total:.12g},"
        f"margin={result.double_cubic_inserted.margin_to_one_third:.12g},"
        "reserve_after_allowance="
        f"{result.reserve_after_declared_allowance:.12g},"
        "status=quarantined_unmasked_coefficient_one_diagnostic"
    )


if __name__ == "__main__":
    main()
