#!/usr/bin/env python3
"""Repair six (1,5,5,1) entries by two physical endpoint factors."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from json import dumps
from math import inf, nextafter, sqrt
from pathlib import Path

from q64_masked_cubic_endpoint_repair import repaired_entries as cubic_repaired_entries
from q64_masked_local_walsh_repair import repaired_entries as local_repaired_entries
from q64_masked_quintic_slice_repair import (
    coefficient_one_dependent_entries,
    exact_quintic_singleton_slice_energies,
    repaired_entries as quintic_repaired_entries,
)


ROOT = Path(__file__).resolve().parents[1]
ORDER = 64
PROFILE = (1, 5, 5, 1)
ProfileSplit = tuple[tuple[int, ...], tuple[int, ...]]


@dataclass(frozen=True)
class MaskedDoubleQuinticEndpointRepair:
    order: int
    candidates: int
    repaired_entries: int
    one_whole_endpoint_entries: int
    oriented_double_one_four_entries: int
    distinct_squared_coefficients: tuple[float, ...]
    minimum_coefficient: float
    maximum_coefficient: float
    remaining_quarantined_entries: int


def previous_repaired_entries() -> frozenset[ProfileSplit]:
    return (
        frozenset(quintic_repaired_entries())
        | frozenset(local_repaired_entries())
        | frozenset(cubic_repaired_entries())
    )


def candidate_entries() -> tuple[ProfileSplit, ...]:
    previous = previous_repaired_entries()
    return tuple(
        sorted(
            entry
            for entry in coefficient_one_dependent_entries()
            if entry not in previous and entry[0] == PROFILE
        )
    )


def endpoint_squared_factor(
    selected: int,
    singleton_on_row: int,
    order: int = ORDER,
) -> Fraction:
    energies = exact_quintic_singleton_slice_energies(order)
    if singleton_on_row not in (0, 1):
        raise ValueError(("singleton split must be binary", singleton_on_row))
    return min(
        order ** (2 * (1 - singleton_on_row)) * energies[selected],
        order ** (2 * singleton_on_row) * energies[5 - selected],
    )


def squared_coefficient(entry: ProfileSplit) -> Fraction:
    profile, split = entry
    if profile != PROFILE:
        raise ValueError(("not a double-quintic endpoint entry", entry))
    return endpoint_squared_factor(split[1], split[0]) * endpoint_squared_factor(
        split[2], split[3]
    )


def repaired_entries() -> tuple[ProfileSplit, ...]:
    return tuple(entry for entry in candidate_entries() if squared_coefficient(entry) <= 1)


def outward_coefficient(entry: ProfileSplit) -> float:
    exact = squared_coefficient(entry)
    value = sqrt(float(exact))
    while Fraction.from_float(value) ** 2 < exact:
        value = nextafter(value, inf)
    return value


def coefficient_map() -> dict[ProfileSplit, float]:
    return {entry: outward_coefficient(entry) for entry in repaired_entries()}


def diagnostic() -> MaskedDoubleQuinticEndpointRepair:
    candidates = candidate_entries()
    repaired = repaired_entries()
    if len(candidates) != 18 or len(repaired) != 6:
        raise AssertionError(("double quintic endpoint inventory", len(candidates), len(repaired)))
    values = tuple(sorted({squared_coefficient(entry) for entry in repaired}))
    coefficients = tuple(outward_coefficient(entry) for entry in repaired)
    return MaskedDoubleQuinticEndpointRepair(
        order=ORDER,
        candidates=len(candidates),
        repaired_entries=len(repaired),
        one_whole_endpoint_entries=sum(
            0 in (entry[1][1], entry[1][2]) or 5 in (entry[1][1], entry[1][2])
            for entry in repaired
        ),
        oriented_double_one_four_entries=sum(
            min(entry[1][1], 5 - entry[1][1]) == 1
            and min(entry[1][2], 5 - entry[1][2]) == 1
            for entry in repaired
        ),
        distinct_squared_coefficients=tuple(map(float, values)),
        minimum_coefficient=min(coefficients),
        maximum_coefficient=max(coefficients),
        remaining_quarantined_entries=(
            len(coefficient_one_dependent_entries())
            - len(previous_repaired_entries())
            - len(repaired)
        ),
    )


def artifact_text(result: MaskedDoubleQuinticEndpointRepair) -> str:
    payload = {
        "schema": "round4_q64_masked_double_quintic_endpoint_repair_v1",
        "result": asdict(result),
        "repaired_registry_entries": [
            {
                "profile": list(profile),
                "split": list(split),
                "squared_coefficient_numerator": squared_coefficient((profile, split)).numerator,
                "squared_coefficient_denominator": squared_coefficient((profile, split)).denominator,
                "outward_coefficient": outward_coefficient((profile, split)),
            }
            for profile, split in repaired_entries()
        ],
        "evidence_label": (
            "arbitrary-correlated-diagonal theorem composing two oriented "
            "complete "
            "physical quintic endpoint row/column factorizations; both "
            "quintic occurrence masks are retained in the exact slices"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.output is not None:
        arguments.output.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 masked double-quintic endpoint repair: "
        f"repaired={result.repaired_entries},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"remaining={result.remaining_quarantined_entries}"
    )


if __name__ == "__main__":
    main()
