#!/usr/bin/env python3
"""Decorated adjacent-row contraction for 16 q64 quintic entries."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from pathlib import Path

from q64_dual_endpoint_schur_insertion import (
    inserted_coefficients as dual_endpoint_inserted_coefficients,
    remaining_quintic_entries as pre_decorated_quintic_entries,
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
from adjacent_balanced_row_slice_contraction import (
    adjacent_balanced_row_coefficient,
)


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class DecoratedAdjacentRowInsertion:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    closed_orbits: int
    degree_twelve_entries: int
    extreme_entries: int
    balanced_entries: int
    fixed_one_cubic_fixed_four_quintic_entries: int
    fixed_two_cubic_fixed_three_quintic_entries: int
    adjacent_row_energy_bound: float
    coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    decorated_row_inserted: OptimizedLedger
    routing_margin_improvement: float
    remaining_quintic_entries: int
    remaining_extreme_entries: int
    remaining_balanced_entries: int
    remaining_quintic_local_proxy: OptimizedLedger
    proxy_reserve_after_declared_allowance: float


def local_complete_row_parameters(
    entry: ProfileSplit,
) -> tuple[int, int, int] | None:
    """Return ``(start, fixed cubic, fixed quintic)`` when the row applies.

    The three consecutive blocks must be singleton--cubic--quintic, in
    either order.  Counts are taken on the side selected by the singleton.
    The original complete-row theorem covers at least one fixed cubic cell
    and three fixed quintic cells; fixing additional cells can only delete
    columns from that row.
    """

    profile, split = entry
    for start in (0, 1):
        outer = 3 if start == 0 else 0
        if profile[outer] != 3 or split[outer] not in (0, 3):
            continue
        degrees = profile[start : start + 3]
        if degrees == (1, 3, 5):
            singleton, cubic, quintic = start, start + 1, start + 2
        elif degrees == (5, 3, 1):
            quintic, cubic, singleton = start, start + 1, start + 2
        else:
            continue
        side = split[singleton]
        if side not in (0, 1):
            raise ValueError(("singleton must be unsplit", entry))
        cubic_fixed = split[cubic] if side else 3 - split[cubic]
        quintic_fixed = split[quintic] if side else 5 - split[quintic]
        if cubic_fixed >= 1 and quintic_fixed >= 3:
            return start, cubic_fixed, quintic_fixed
    return None


def decorated_adjacent_entries() -> tuple[ProfileSplit, ...]:
    return tuple(
        entry
        for entry in pre_decorated_quintic_entries()
        if local_complete_row_parameters(entry) is not None
    )


def remaining_quintic_entries() -> tuple[ProfileSplit, ...]:
    closed = set(decorated_adjacent_entries())
    return tuple(
        entry
        for entry in pre_decorated_quintic_entries()
        if entry not in closed
    )


def decorated_adjacent_coefficient() -> float:
    return adjacent_balanced_row_coefficient(DIMENSION)


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = dual_endpoint_inserted_coefficients()
    coefficient = decorated_adjacent_coefficient()
    for entry in decorated_adjacent_entries():
        result[entry] = coefficient
    return result


def remaining_quintic_local_proxy_coefficients() -> dict[ProfileSplit, float]:
    # Import lazily so the earlier theorem remains independent of this one.
    from q64_dual_endpoint_schur_insertion import local_slice_coefficients

    result = inserted_coefficients()
    extreme, balanced = local_slice_coefficients()
    for entry in remaining_quintic_entries():
        result[entry] = (
            extreme if quintic_split_depth(entry) == 1 else balanced
        )
    return result


def diagnostic() -> DecoratedAdjacentRowInsertion:
    entries = decorated_adjacent_entries()
    remaining = remaining_quintic_entries()
    coefficient = decorated_adjacent_coefficient()
    parameters = tuple(local_complete_row_parameters(entry) for entry in entries)
    previous = optimize(
        mapped_coefficients=dual_endpoint_inserted_coefficients()
    )
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    proxy = optimize(
        mapped_coefficients=remaining_quintic_local_proxy_coefficients()
    )
    previous_proved = 276
    return DecoratedAdjacentRowInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        closed_orbits=len(entries) // 4,
        degree_twelve_entries=sum(sum(entry[0]) == 12 for entry in entries),
        extreme_entries=sum(quintic_split_depth(entry) == 1 for entry in entries),
        balanced_entries=sum(quintic_split_depth(entry) == 2 for entry in entries),
        fixed_one_cubic_fixed_four_quintic_entries=sum(
            parameter is not None and parameter[1:] == (1, 4)
            for parameter in parameters
        ),
        fixed_two_cubic_fixed_three_quintic_entries=sum(
            parameter is not None and parameter[1:] == (2, 3)
            for parameter in parameters
        ),
        adjacent_row_energy_bound=coefficient**2,
        coefficient=coefficient,
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        decorated_row_inserted=inserted,
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


def artifact_text(result: DecoratedAdjacentRowInsertion) -> str:
    payload = {
        "schema": "round4_q64_decorated_adjacent_row_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal complete-row contraction inherited from the "
            "proved adjacent cubic-quintic row; the unsplit fourth block is "
            "a unit cross-Gram Schur multiplier; floating q64 Perron "
            "insertion; one batch only"
        ),
        "remaining_open": (
            "596 balanced entries, including 136 split-cubic/split-quintic "
            "entries, plus interval certification and the adaptive lift"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic decorated-row insertion",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_decorated_adjacent_row_insertion.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 decorated adjacent-row insertion: "
        f"entries={result.closed_entries},"
        f"orbits={result.closed_orbits},"
        f"coefficient={result.coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.decorated_row_inserted.total:.12g},"
        f"margin={result.decorated_row_inserted.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        f"proxy_total={result.remaining_quintic_local_proxy.total:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
