#!/usr/bin/env python3
"""Dual endpoint-slice Schur contraction for 12 q64 quintic entries."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from json import dumps
from math import inf, nextafter, sqrt
from pathlib import Path

from q64_fixed_pair_adjacent_row_contraction import (
    inserted_coefficients as fixed_pair_inserted_coefficients,
)
from q64_fixed_singleton_pair_contraction import fixed_singleton_pair_entries
from q64_paper_target_gate import (
    DIMENSION,
    MODES,
    ORDER,
    RESERVE_TARGET,
    OptimizedLedger,
    ProfileSplit,
    optimize,
)
from q64_post_universal_quintic_gate import (
    quintic_entries,
    quintic_split_depth,
)
from q64_reversed_middle_pair_contraction import reversed_middle_pair_entries
from q64_shifted_middle_pair_contraction import shifted_middle_pair_entries
import occupation_compatible_sector_optimization as occupation


ROOT = Path(__file__).resolve().parents[1]


def cubic_fixed_pair_energy(order: int = ORDER) -> Fraction:
    """Exact physical cubic endpoint energy through a fixed pair."""

    q = order
    return Fraction(q * q - 2 * q + 2, q * q * (q - 1))


def quintic_fixed_triple_energy(order: int = ORDER) -> Fraction:
    """Exact physical quintic endpoint energy through a fixed triple."""

    q = order
    w0 = Fraction(1, q * q)
    w1 = Fraction(1, q * q * (q - 1) ** 2)
    return (
        (q - 3) * (q - 4) // 2 * w0
        + q * (q - 1) * (w0 + (q - 4) * w1)
        + (q - 1) * q * (q - 1) // 2 * w1
    )


def outward_sqrt(value: Fraction) -> float:
    """Return a binary64 upper bound for the square root of ``value``."""

    result = sqrt(value.numerator / value.denominator)
    while Fraction.from_float(result) ** 2 < value:
        result = nextafter(result, inf)
    return result


@dataclass(frozen=True)
class DualEndpointSchurInsertion:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    closed_orbits: int
    degree_ten_entries: int
    degree_twelve_entries: int
    extreme_entries: int
    balanced_entries: int
    cubic_fixed_pair_energy: float
    cubic_schur_factor: float
    balanced_quintic_fixed_triple_energy: float
    balanced_quintic_schur_factor: float
    balanced_coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    dual_endpoint_inserted: OptimizedLedger
    routing_margin_spent: float
    remaining_quintic_entries: int
    remaining_extreme_entries: int
    remaining_balanced_entries: int
    remaining_quintic_local_proxy: OptimizedLedger
    proxy_reserve_after_declared_allowance: float


def pre_dual_quintic_entries() -> tuple[ProfileSplit, ...]:
    closed = (
        set(shifted_middle_pair_entries())
        | set(reversed_middle_pair_entries())
        | set(fixed_singleton_pair_entries())
        | set(fixed_pair_adjacent_entries())
    )
    return tuple(entry for entry in quintic_entries() if entry not in closed)


def fixed_pair_adjacent_entries() -> tuple[ProfileSplit, ...]:
    # Local import avoids making the earlier theorem depend on this module.
    from q64_fixed_pair_adjacent_row_contraction import (
        fixed_pair_adjacent_entries as entries,
    )

    return entries()


def split_cubic_index(entry: ProfileSplit) -> int:
    profile, split = entry
    indices = tuple(
        index
        for index, (degree, selected) in enumerate(
            zip(profile, split, strict=True)
        )
        if degree == 3 and selected in (1, 2)
    )
    if len(indices) != 1:
        raise ValueError(("expected one split cubic", entry, indices))
    return indices[0]


def has_favorable_cubic_singleton(entry: ProfileSplit) -> bool:
    profile, split = entry
    index = split_cubic_index(entry)
    pair_side = split[index] - 1
    return any(
        profile[neighbor] == 1 and split[neighbor] == pair_side
        for neighbor in (index - 1, index + 1)
        if 0 <= neighbor < 4
    )


def has_favorable_quintic_singleton(entry: ProfileSplit) -> bool:
    profile, split = entry
    index = profile.index(5)
    majority_side = int(split[index] > 2)
    return any(
        profile[neighbor] == 1 and split[neighbor] == majority_side
        for neighbor in (index - 1, index + 1)
        if 0 <= neighbor < 4
    )


def dual_endpoint_entries() -> tuple[ProfileSplit, ...]:
    return tuple(
        entry
        for entry in pre_dual_quintic_entries()
        if has_favorable_cubic_singleton(entry)
        and has_favorable_quintic_singleton(entry)
    )


def remaining_quintic_entries() -> tuple[ProfileSplit, ...]:
    closed = set(dual_endpoint_entries())
    return tuple(
        entry for entry in pre_dual_quintic_entries() if entry not in closed
    )


def local_slice_coefficients() -> tuple[float, float]:
    cubic = cubic_fixed_pair_energy()
    quintic = quintic_fixed_triple_energy()
    extreme = Fraction(ORDER * ORDER - 4, ORDER * ORDER)
    return outward_sqrt(cubic * extreme), outward_sqrt(cubic * quintic)


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = fixed_pair_inserted_coefficients()
    extreme, balanced = local_slice_coefficients()
    for entry in dual_endpoint_entries():
        result[entry] = (
            extreme if quintic_split_depth(entry) == 1 else balanced
        )
    return result


def remaining_quintic_local_proxy_coefficients() -> dict[ProfileSplit, float]:
    result = inserted_coefficients()
    extreme, balanced = local_slice_coefficients()
    for entry in remaining_quintic_entries():
        result[entry] = (
            extreme if quintic_split_depth(entry) == 1 else balanced
        )
    return result


def diagnostic() -> DualEndpointSchurInsertion:
    entries = dual_endpoint_entries()
    remaining = remaining_quintic_entries()
    cubic = occupation.endpoint_singleton_slice_energies(ORDER)[2]
    quintic = occupation.endpoint_quintic_singleton_slice_energies(ORDER)
    _, balanced = local_slice_coefficients()
    previous = optimize(mapped_coefficients=fixed_pair_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    proxy = optimize(
        mapped_coefficients=remaining_quintic_local_proxy_coefficients()
    )
    previous_proved = 264
    return DualEndpointSchurInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        closed_orbits=len(entries) // 4,
        degree_ten_entries=sum(sum(entry[0]) == 10 for entry in entries),
        degree_twelve_entries=sum(sum(entry[0]) == 12 for entry in entries),
        extreme_entries=sum(quintic_split_depth(entry) == 1 for entry in entries),
        balanced_entries=sum(quintic_split_depth(entry) == 2 for entry in entries),
        cubic_fixed_pair_energy=cubic,
        cubic_schur_factor=sqrt(cubic),
        balanced_quintic_fixed_triple_energy=quintic[3],
        balanced_quintic_schur_factor=sqrt(quintic[3]),
        balanced_coefficient=balanced,
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        dual_endpoint_inserted=inserted,
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


def artifact_text(result: DualEndpointSchurInsertion) -> str:
    cubic = cubic_fixed_pair_energy()
    quintic = quintic_fixed_triple_energy()
    squared = cubic * quintic
    payload = {
        "schema": "round4_q64_dual_endpoint_schur_insertion_v2",
        "result": asdict(result),
        "exact_certificate": {
            "cubic_energy_numerator": cubic.numerator,
            "cubic_energy_denominator": cubic.denominator,
            "quintic_energy_numerator": quintic.numerator,
            "quintic_energy_denominator": quintic.denominator,
            "squared_coefficient_numerator": squared.numerator,
            "squared_coefficient_denominator": squared.denominator,
            "outward_coefficient": local_slice_coefficients()[1],
        },
        "registry_entries": [
            {"profile": list(profile), "split": list(split)}
            for profile, split in dual_endpoint_entries()
        ],
        "evidence_label": (
            "arbitrary-diagonal trace-class Schur-multiplier composition "
            "using exact cubic fixed-pair and quintic fixed-triple endpoint "
            "slices, independently checked completed-link Gram lifts, and "
            "an outward-rounded exact-rational coefficient; historical "
            "floating q64 insertion retained only for provenance; one batch only"
        ),
        "historical_checkpoint": (
            "At the original insertion this left 612 balanced entries, including "
            "152 split-cubic/split-quintic entries. Subsequent theorems close "
            "that inventory; the current registry is 888 of 888 one-batch "
            "entries, and only the adaptive lift remains open."
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic dual-endpoint Schur insertion",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_dual_endpoint_schur_insertion.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 dual-endpoint Schur insertion: "
        f"entries={result.closed_entries},"
        f"orbits={result.closed_orbits},"
        f"coefficient={result.balanced_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.dual_endpoint_inserted.total:.12g},"
        f"margin={result.dual_endpoint_inserted.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        f"proxy_total={result.remaining_quintic_local_proxy.total:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
