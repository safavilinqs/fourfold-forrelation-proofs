#!/usr/bin/env python3
"""Close the 96 low-impact septimic q64 entries by the global Gram bound.

For any probability law on sign vectors and any fixed occurrence split, the
moment matrix is a cross Gram matrix of unit-modulus character features.
Arbitrary diagonal weighting therefore has nuclear norm at most the square
root of the two total masses.  The resulting universal coefficient one is
coarse, but the 96 septimic entries in the leading remaining class have low
enough Perron impact that inserting it retains the declared q64 reserve.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from pathlib import Path

from q64_block_coherent_contraction import block_coherent_entries
from q64_chain_aware_insertion import (
    chain_aware_entries,
    inserted_coefficients as chain_inserted_coefficients,
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
from q64_remaining_class_gates import contraction_class, remaining_entries


ROOT = Path(__file__).resolve().parents[1]
UNIVERSAL_GRAM_COEFFICIENT = 1.0


@dataclass(frozen=True)
class UniversalSeptimicInsertion:
    order: int
    dimension: int
    sign_modes: int
    universal_coefficient: float
    septimic_entries: int
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    septimic_inserted: OptimizedLedger
    margin_spent: float
    reserve_after_declared_allowance: float


def septimic_entries() -> tuple[ProfileSplit, ...]:
    return tuple(
        entry
        for entry in remaining_entries()
        if contraction_class(entry) == "one_split_cubic_one_split_higher"
        and 7 in entry[0]
    )


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = chain_inserted_coefficients()
    for entry in septimic_entries():
        result[entry] = UNIVERSAL_GRAM_COEFFICIENT
    return result


def diagnostic() -> UniversalSeptimicInsertion:
    entries = septimic_entries()
    previous_proved = len(block_coherent_entries()) + len(chain_aware_entries())
    previous = optimize(mapped_coefficients=chain_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    return UniversalSeptimicInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        universal_coefficient=UNIVERSAL_GRAM_COEFFICIENT,
        septimic_entries=len(entries),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        septimic_inserted=inserted,
        margin_spent=(
            previous.margin_to_one_third - inserted.margin_to_one_third
        ),
        reserve_after_declared_allowance=(
            inserted.margin_to_one_third - RESERVE_TARGET
        ),
    )


def artifact_text(result: UniversalSeptimicInsertion) -> str:
    payload = {
        "schema": "round4_q64_universal_septimic_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "quarantined unmasked cross-Gram diagnostic for 96 septimic "
            "entries; it omits cross-cut distinctness masks and is not a theorem"
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
        help="write the deterministic universal septimic insertion",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_universal_septimic_insertion.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 universal septimic insertion: "
        f"entries={result.septimic_entries},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.septimic_inserted.total:.12g},"
        f"margin={result.septimic_inserted.margin_to_one_third:.12g},"
        "reserve_after_allowance="
        f"{result.reserve_after_declared_allowance:.12g},"
        "status=quarantined_unmasked_coefficient_one_diagnostic"
    )


if __name__ == "__main__":
    main()
