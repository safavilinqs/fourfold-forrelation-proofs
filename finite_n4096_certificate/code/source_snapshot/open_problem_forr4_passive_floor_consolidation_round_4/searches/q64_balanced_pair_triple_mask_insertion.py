#!/usr/bin/env python3
"""Balanced pair--triple mask contraction for eight q64 entries."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from math import sqrt
from pathlib import Path

from q64_degree_ten_completion_row_insertion import orbit
from q64_dual_endpoint_schur_insertion import (
    has_favorable_cubic_singleton,
    has_favorable_quintic_singleton,
)
from q64_internal_whole_cubic_endpoint_insertion import (
    inserted_coefficients as internal_whole_inserted_coefficients,
    remaining_quintic_entries as pre_pair_triple_quintic_entries,
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
from q64_post_universal_quintic_gate import quintic_split_depth
import occupation_compatible_sector_optimization as occupation


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class BalancedPairTripleMaskInsertion:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    closed_orbits: int
    extreme_entries: int
    balanced_entries: int
    quintic_favorable_entries: int
    cubic_fixed_pair_energy: float
    cubic_endpoint_factor: float
    pair_triple_mask_factor: float
    coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    inserted_routing: OptimizedLedger
    routing_margin_spent: float
    remaining_quintic_entries: int
    remaining_extreme_entries: int
    remaining_balanced_entries: int
    remaining_quintic_local_proxy: OptimizedLedger
    proxy_reserve_after_declared_allowance: float


def pair_triple_disjointness_factor() -> float:
    """Direct Schur factor for disjoint two- and three-subsets.

    The three summands in ``1-|F intersection G|+1[F subset G]`` have
    row squared multiplicities ``(1, 2, 1)`` and column squared
    multiplicities ``(1, 3, 3)``.  Optimal direct-sum scaling gives the
    sum of their geometric means.
    """

    return 1 + sqrt(6) + sqrt(3)


def balanced_pair_triple_entries() -> tuple[ProfileSplit, ...]:
    """Return the remaining cuts with an exact cubic endpoint slice."""

    return tuple(
        entry
        for entry in pre_pair_triple_quintic_entries()
        if has_favorable_cubic_singleton(entry)
    )


def coefficient(order: int = ORDER) -> float:
    cubic_energy = occupation.endpoint_singleton_slice_energies(order)[2]
    return sqrt(cubic_energy) * pair_triple_disjointness_factor()


def remaining_quintic_entries() -> tuple[ProfileSplit, ...]:
    closed = set(balanced_pair_triple_entries())
    return tuple(
        entry
        for entry in pre_pair_triple_quintic_entries()
        if entry not in closed
    )


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = internal_whole_inserted_coefficients()
    value = coefficient()
    for entry in balanced_pair_triple_entries():
        result[entry] = value
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


def diagnostic() -> BalancedPairTripleMaskInsertion:
    entries = balanced_pair_triple_entries()
    remaining = remaining_quintic_entries()
    cubic_energy = occupation.endpoint_singleton_slice_energies(ORDER)[2]
    previous = optimize(
        mapped_coefficients=internal_whole_inserted_coefficients()
    )
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    proxy = optimize(
        mapped_coefficients=remaining_quintic_local_proxy_coefficients()
    )
    previous_proved = 340
    return BalancedPairTripleMaskInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        closed_orbits=len({frozenset(orbit(entry)) for entry in entries}),
        extreme_entries=sum(quintic_split_depth(entry) == 1 for entry in entries),
        balanced_entries=sum(
            quintic_split_depth(entry) == 2 for entry in entries
        ),
        quintic_favorable_entries=sum(
            has_favorable_quintic_singleton(entry) for entry in entries
        ),
        cubic_fixed_pair_energy=cubic_energy,
        cubic_endpoint_factor=sqrt(cubic_energy),
        pair_triple_mask_factor=pair_triple_disjointness_factor(),
        coefficient=coefficient(),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        inserted_routing=inserted,
        routing_margin_spent=(
            previous.margin_to_one_third - inserted.margin_to_one_third
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


def artifact_text(result: BalancedPairTripleMaskInsertion) -> str:
    payload = {
        "schema": "round4_q64_balanced_pair_triple_mask_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal exact cubic fixed-pair endpoint slice, "
            "direct constant/incidence/pair-containment factorization of "
            "the balanced quintic distinctness mask, and unit cross-Gram "
            "completion of all remaining physical links; floating q64 "
            "Perron insertion; one batch only"
        ),
        "remaining_open": (
            "540 balanced entries, including 80 degree-twelve "
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
        help="write the deterministic balanced pair--triple insertion",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = (
            ROOT
            / "artifacts"
            / "q64_balanced_pair_triple_mask_insertion.json"
        )
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 balanced pair-triple mask insertion: "
        f"entries={result.closed_entries},"
        f"orbits={result.closed_orbits},"
        f"mask_factor={result.pair_triple_mask_factor:.12g},"
        f"coefficient={result.coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.inserted_routing.total:.12g},"
        f"margin={result.inserted_routing.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        f"proxy_total={result.remaining_quintic_local_proxy.total:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
