#!/usr/bin/env python3
"""Internal whole-cubic endpoint contraction for 16 q64 entries."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from pathlib import Path

from q64_degree_ten_completion_row_insertion import orbit
from q64_last_degree_ten_chain_insertion import (
    inserted_coefficients as last_degree_ten_inserted_coefficients,
    remaining_quintic_entries as pre_internal_whole_quintic_entries,
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
from leading_balanced_disjointness_contraction import (
    disjointness_schur_factor,
)


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class InternalWholeCubicEndpointInsertion:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    closed_orbits: int
    extreme_entries: int
    balanced_entries: int
    cubic_distinctness_factor: float
    extreme_quintic_distinctness_factor: float
    balanced_quintic_distinctness_factor: float
    endpoint_factor: float
    extreme_coefficient: float
    balanced_coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    endpoint_inserted: OptimizedLedger
    routing_change: float
    remaining_quintic_entries: int
    remaining_extreme_entries: int
    remaining_balanced_entries: int
    remaining_quintic_local_proxy: OptimizedLedger
    proxy_reserve_after_declared_allowance: float


def same_side_adjacent_whole_cubic(entry: ProfileSplit) -> bool:
    profile, split = entry
    singleton = profile.index(1)
    whole_cubics = tuple(
        index
        for index, (degree, selected) in enumerate(
            zip(profile, split, strict=True)
        )
        if degree == 3 and selected in (0, 3)
    )
    if len(whole_cubics) != 1:
        return False
    whole = whole_cubics[0]
    return abs(whole - singleton) == 1 and split[whole] == 3 * split[singleton]


def internal_whole_cubic_entries() -> tuple[ProfileSplit, ...]:
    return tuple(
        entry
        for entry in pre_internal_whole_quintic_entries()
        if same_side_adjacent_whole_cubic(entry)
    )


def distinctness_factors(
    dimension: int = DIMENSION,
) -> tuple[float, float, float]:
    cubic = disjointness_schur_factor(dimension, 2)
    extreme = disjointness_schur_factor(dimension, 4)
    balanced = disjointness_schur_factor(dimension, 3) ** 2
    return cubic, extreme, balanced


def coefficients(order: int = ORDER) -> tuple[float, float]:
    cubic, extreme, balanced = distinctness_factors(order * order)
    return cubic * extreme / order, cubic * balanced / order


def remaining_quintic_entries() -> tuple[ProfileSplit, ...]:
    closed = set(internal_whole_cubic_entries())
    return tuple(
        entry
        for entry in pre_internal_whole_quintic_entries()
        if entry not in closed
    )


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = last_degree_ten_inserted_coefficients()
    extreme, balanced = coefficients()
    for entry in internal_whole_cubic_entries():
        result[entry] = (
            extreme if quintic_split_depth(entry) == 1 else balanced
        )
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


def diagnostic() -> InternalWholeCubicEndpointInsertion:
    entries = internal_whole_cubic_entries()
    remaining = remaining_quintic_entries()
    cubic, extreme_factor, balanced_factor = distinctness_factors()
    extreme_coefficient, balanced_coefficient = coefficients()
    previous = optimize(
        mapped_coefficients=last_degree_ten_inserted_coefficients()
    )
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    proxy = optimize(
        mapped_coefficients=remaining_quintic_local_proxy_coefficients()
    )
    previous_proved = 324
    return InternalWholeCubicEndpointInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        closed_orbits=len({frozenset(orbit(entry)) for entry in entries}),
        extreme_entries=sum(quintic_split_depth(entry) == 1 for entry in entries),
        balanced_entries=sum(quintic_split_depth(entry) == 2 for entry in entries),
        cubic_distinctness_factor=cubic,
        extreme_quintic_distinctness_factor=extreme_factor,
        balanced_quintic_distinctness_factor=balanced_factor,
        endpoint_factor=1 / ORDER,
        extreme_coefficient=extreme_coefficient,
        balanced_coefficient=balanced_coefficient,
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        endpoint_inserted=inserted,
        routing_change=inserted.total - previous.total,
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


def artifact_text(result: InternalWholeCubicEndpointInsertion) -> str:
    payload = {
        "schema": "round4_q64_internal_whole_cubic_endpoint_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal product of completed physical link kernels, "
            "explicit cubic/quintic distinctness Schur factors, and a "
            "same-side whole-cubic endpoint moment bounded by 1/q; floating "
            "q64 Perron insertion; one batch only"
        ),
        "remaining_open": (
            "548 balanced entries, including 88 degree-twelve "
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
        help="write the deterministic internal whole-cubic endpoint insertion",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = (
            ROOT
            / "artifacts"
            / "q64_internal_whole_cubic_endpoint_insertion.json"
        )
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 internal whole-cubic endpoint insertion: "
        f"entries={result.closed_entries},"
        f"orbits={result.closed_orbits},"
        f"coefficients={result.extreme_coefficient:.12g}/"
        f"{result.balanced_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.endpoint_inserted.total:.12g},"
        f"margin={result.endpoint_inserted.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        f"proxy_total={result.remaining_quintic_local_proxy.total:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
