#!/usr/bin/env python3
"""Final degree-ten chain contraction for four q64 quintic entries."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from math import sqrt
from pathlib import Path

from q64_degree_ten_completion_row_insertion import orbit
from q64_paper_target_gate import (
    DIMENSION,
    MODES,
    ORDER,
    RESERVE_TARGET,
    OptimizedLedger,
    ProfileSplit,
    optimize,
)
from q64_post_universal_quintic_gate import quintic_split_depth
from q64_whole_cubic_decorated_row_insertion import (
    inserted_coefficients as whole_cubic_inserted_coefficients,
    remaining_quintic_entries as pre_last_degree_ten_quintic_entries,
)
from leading_balanced_disjointness_contraction import (
    disjointness_schur_factor,
)


ROOT = Path(__file__).resolve().parents[1]
TARGET = ((3, 1, 1, 5), (1, 1, 0, 3))


@dataclass(frozen=True)
class LastDegreeTenChainInsertion:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    closed_orbits: int
    cubic_distinctness_factor: float
    quintic_completion_factor: float
    walsh_chain_factor: float
    coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    last_degree_ten_inserted: OptimizedLedger
    routing_margin_improvement: float
    remaining_quintic_entries: int
    remaining_extreme_entries: int
    remaining_balanced_entries: int
    remaining_quintic_local_proxy: OptimizedLedger
    proxy_reserve_after_declared_allowance: float


def last_degree_ten_entries() -> tuple[ProfileSplit, ...]:
    return orbit(TARGET)


def last_degree_ten_coefficient(
    order: int = ORDER,
) -> float:
    dimension = order * order
    cubic = disjointness_schur_factor(dimension, 2)
    quintic = 1 + sqrt(2)
    return cubic * quintic / order


def remaining_quintic_entries() -> tuple[ProfileSplit, ...]:
    closed = set(last_degree_ten_entries())
    return tuple(
        entry
        for entry in pre_last_degree_ten_quintic_entries()
        if entry not in closed
    )


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = whole_cubic_inserted_coefficients()
    coefficient = last_degree_ten_coefficient()
    for entry in last_degree_ten_entries():
        result[entry] = coefficient
    return result


def remaining_quintic_local_proxy_coefficients() -> dict[ProfileSplit, float]:
    from q64_dual_endpoint_schur_insertion import local_slice_coefficients

    result = inserted_coefficients()
    extreme, balanced = local_slice_coefficients()
    for entry in remaining_quintic_entries():
        result[entry] = (
            extreme if quintic_split_depth(entry) == 1 else balanced
        )
    return result


def diagnostic() -> LastDegreeTenChainInsertion:
    entries = last_degree_ten_entries()
    remaining = remaining_quintic_entries()
    previous = optimize(mapped_coefficients=whole_cubic_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    proxy = optimize(
        mapped_coefficients=remaining_quintic_local_proxy_coefficients()
    )
    previous_proved = 320
    return LastDegreeTenChainInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        closed_orbits=len(entries) // 4,
        cubic_distinctness_factor=disjointness_schur_factor(DIMENSION, 2),
        quintic_completion_factor=1 + sqrt(2),
        walsh_chain_factor=1 / ORDER,
        coefficient=last_degree_ten_coefficient(),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        last_degree_ten_inserted=inserted,
        routing_margin_improvement=(
            inserted.margin_to_one_third - previous.margin_to_one_third
        ),
        remaining_quintic_entries=len(remaining),
        remaining_extreme_entries=sum(
            quintic_split_depth(entry) == 1 for entry in remaining
        ),
        remaining_balanced_entries=sum(
            quintic_split_depth(entry) == 2 for entry in remaining
        ),
        remaining_quintic_local_proxy=proxy,
        proxy_reserve_after_declared_allowance=(
            proxy.margin_to_one_third - RESERVE_TARGET
        ),
    )


def artifact_text(result: LastDegreeTenChainInsertion) -> str:
    payload = {
        "schema": "round4_q64_last_degree_ten_chain_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal whole-chain Schur contraction using a "
            "completed cubic endpoint plus a fixed-singleton/pair "
            "distinctness factor, the transposed completed quintic endpoint "
            "bound, and an exact 1/q Walsh-chain collapse; floating q64 "
            "Perron insertion; one batch only"
        ),
        "remaining_open": (
            "564 balanced entries, including 104 degree-twelve "
            "split-cubic/split-quintic entries, plus interval certification "
            "and the adaptive lift"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic final degree-ten chain insertion",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_last_degree_ten_chain_insertion.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 final degree-ten chain insertion: "
        f"entries={result.closed_entries},"
        f"coefficient={result.coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.last_degree_ten_inserted.total:.12g},"
        f"margin={result.last_degree_ten_inserted.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        f"proxy_total={result.remaining_quintic_local_proxy.total:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
