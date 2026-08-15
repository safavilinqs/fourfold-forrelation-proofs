#!/usr/bin/env python3
"""Insert the ten inherited chain-aware orbit theorems at q=64.

The Round 3 theorem formulas are dimension-parameterized even though the
original finite-size ledger used q=32.  At N=4096 they close forty internally
split entries, disjoint from the seventy block-coherent entries already
closed in Round 4.  This module evaluates those exact formulas, checks the
partition, and reoptimizes the q64 occupation ledger.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent / "open_problem_forr4_passive_floor_consolidation_round_3" / "searches"
)
sys.path.insert(0, str(ROUND3_SEARCHES))

from finite_size_passive_ledger import theorem_families  # noqa: E402
from q64_block_coherent_contraction import (  # noqa: E402
    block_coherent_entries,
    inserted_coefficients as block_inserted_coefficients,
)
from q64_paper_target_gate import (  # noqa: E402
    DIMENSION,
    MODES,
    ORDER,
    OptimizedLedger,
    optimize,
)


Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class FamilyRow:
    name: str
    coefficient: float
    entries: int
    source: str


@dataclass(frozen=True)
class ChainAwareInsertion:
    order: int
    dimension: int
    sign_modes: int
    theorem_families: int
    chain_aware_entries: int
    overlap_with_block_coherent: int
    total_proved_open_entries: int
    remaining_open_entries: int
    minimum_coefficient: float
    minimum_family: str
    maximum_coefficient: float
    maximum_family: str
    families: tuple[FamilyRow, ...]
    block_coherent_only: OptimizedLedger
    chain_aware_inserted: OptimizedLedger
    margin_improvement: float


def chain_aware_entries() -> dict[ProfileSplit, float]:
    result: dict[ProfileSplit, float] = {}
    for family in theorem_families(DIMENSION):
        for entry in family.entries:
            if entry in result:
                raise AssertionError(("chain-aware theorem overlap", entry))
            result[entry] = family.coefficient
    return result


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = block_inserted_coefficients()
    result.update(chain_aware_entries())
    return result


def diagnostic() -> ChainAwareInsertion:
    families = theorem_families(DIMENSION)
    entries = chain_aware_entries()
    block_entries = set(block_coherent_entries())
    overlap = block_entries.intersection(entries)
    if overlap:
        raise AssertionError(("block/chain theorem overlap", overlap))
    minimum = min(families, key=lambda family: family.coefficient)
    maximum = max(families, key=lambda family: family.coefficient)
    previous = optimize(mapped_coefficients=block_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    return ChainAwareInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        theorem_families=len(families),
        chain_aware_entries=len(entries),
        overlap_with_block_coherent=len(overlap),
        total_proved_open_entries=len(block_entries) + len(entries),
        remaining_open_entries=888 - len(block_entries) - len(entries),
        minimum_coefficient=minimum.coefficient,
        minimum_family=minimum.name,
        maximum_coefficient=maximum.coefficient,
        maximum_family=maximum.name,
        families=tuple(
            FamilyRow(
                name=family.name,
                coefficient=family.coefficient,
                entries=len(family.entries),
                source=family.source,
            )
            for family in families
        ),
        block_coherent_only=previous,
        chain_aware_inserted=inserted,
        margin_improvement=(
            inserted.margin_to_one_third - previous.margin_to_one_third
        ),
    )


def artifact_text(result: ChainAwareInsertion) -> str:
    payload = {
        "schema": "round4_q64_chain_aware_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "ten inherited arbitrary-diagonal orbit theorems evaluated at "
            "q=64 and combined with exact block-coherent coefficients; "
            "floating Perron insertion; one batch only"
        ),
        "remaining_open": (
            "778 internally split balanced entries plus interval certification "
            "and the adaptive lift"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic q64 chain-aware insertion",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_chain_aware_insertion.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 chain-aware insertion: "
        f"families={result.theorem_families},"
        f"new_entries={result.chain_aware_entries},"
        f"proved_entries={result.total_proved_open_entries},"
        f"remaining={result.remaining_open_entries},"
        f"coefficient_min={result.minimum_coefficient:.12g},"
        f"coefficient_max={result.maximum_coefficient:.12g},"
        f"total={result.chain_aware_inserted.total:.12g},"
        f"margin={result.chain_aware_inserted.margin_to_one_third:.12g},"
        f"margin_gain={result.margin_improvement:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
