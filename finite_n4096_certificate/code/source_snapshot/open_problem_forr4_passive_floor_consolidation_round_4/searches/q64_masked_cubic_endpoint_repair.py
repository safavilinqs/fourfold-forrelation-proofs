#!/usr/bin/env python3
"""Repair 12 masked entries by a physical cubic endpoint slice."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from json import dumps
from math import inf, nextafter, sqrt
from pathlib import Path

from q64_masked_local_walsh_repair import repaired_entries as local_repaired_entries
from q64_masked_quintic_slice_repair import (
    coefficient_one_dependent_entries,
    repaired_entries as quintic_repaired_entries,
)


ROOT = Path(__file__).resolve().parents[1]
ORDER = 64
ProfileSplit = tuple[tuple[int, ...], tuple[int, ...]]


@dataclass(frozen=True)
class MaskedCubicEndpointRepair:
    order: int
    candidate_entries: int
    repaired_entries: int
    cubic_septimic_entries: int
    recovered_cubic_quintic_entries: int
    exact_cubic_endpoint_squared_factor: float
    residual_mask_factor: float
    squared_coefficient: float
    coefficient: float
    remaining_quarantined_entries: int


def cubic_endpoint_squared_factor(order: int = ORDER) -> Fraction:
    q = order
    return Fraction(q * q - 2 * q + 2, q * q * (q - 1))


def previous_repaired_entries() -> frozenset[ProfileSplit]:
    return frozenset(quintic_repaired_entries()) | frozenset(local_repaired_entries())


def entry_kind(entry: ProfileSplit) -> str | None:
    profile, split = entry
    endpoint_cubics = tuple(
        cubic
        for singleton, cubic in ((0, 1), (3, 2))
        if profile[singleton] == 1
        and profile[cubic] == 3
        and split[cubic] in (1, 2)
    )
    if not endpoint_cubics:
        return None
    if sorted(profile) == [1, 1, 3, 7]:
        septimic = profile.index(7)
        if min(split[septimic], 7 - split[septimic]) == 2:
            return "cubic_septimic"
    if sorted(profile) == [1, 3, 3, 5]:
        quintic = profile.index(5)
        if min(split[quintic], 5 - split[quintic]) == 1:
            return "recovered_cubic_quintic"
    return None


def candidate_entries() -> tuple[ProfileSplit, ...]:
    previous = previous_repaired_entries()
    return tuple(
        sorted(
            entry
            for entry in coefficient_one_dependent_entries()
            if entry not in previous and entry_kind(entry) is not None
        )
    )


def squared_coefficient(order: int = ORDER) -> Fraction:
    # gamma_12 < 5/2 and gamma_14 = 3, while gamma_25 < 15/2.
    # Both structural cases therefore leave residual mask factor at most 15/2.
    return Fraction(225, 4) * cubic_endpoint_squared_factor(order)


def coefficient(order: int = ORDER) -> float:
    exact = squared_coefficient(order)
    value = sqrt(float(exact))
    while Fraction.from_float(value) ** 2 < exact:
        value = nextafter(value, inf)
    return value


def repaired_entries() -> tuple[ProfileSplit, ...]:
    entries = candidate_entries()
    return entries if squared_coefficient() < 1 else ()


def coefficient_map() -> dict[ProfileSplit, float]:
    return {entry: coefficient() for entry in repaired_entries()}


def diagnostic() -> MaskedCubicEndpointRepair:
    entries = repaired_entries()
    if len(entries) != 12:
        raise AssertionError(("masked cubic endpoint inventory", len(entries)))
    kinds = tuple(entry_kind(entry) for entry in entries)
    return MaskedCubicEndpointRepair(
        order=ORDER,
        candidate_entries=len(candidate_entries()),
        repaired_entries=len(entries),
        cubic_septimic_entries=kinds.count("cubic_septimic"),
        recovered_cubic_quintic_entries=kinds.count("recovered_cubic_quintic"),
        exact_cubic_endpoint_squared_factor=float(cubic_endpoint_squared_factor()),
        residual_mask_factor=7.5,
        squared_coefficient=float(squared_coefficient()),
        coefficient=coefficient(),
        remaining_quarantined_entries=(
            len(coefficient_one_dependent_entries())
            - len(previous_repaired_entries())
            - len(entries)
        ),
    )


def artifact_text(result: MaskedCubicEndpointRepair) -> str:
    payload = {
        "schema": "round4_q64_masked_cubic_endpoint_repair_v1",
        "result": asdict(result),
        "repaired_registry_entries": [
            {
                "profile": list(profile),
                "split": list(split),
                "kind": entry_kind((profile, split)),
                "squared_coefficient_numerator": squared_coefficient().numerator,
                "squared_coefficient_denominator": squared_coefficient().denominator,
                "outward_coefficient": coefficient(),
            }
            for profile, split in repaired_entries()
        ],
        "evidence_label": (
            "arbitrary-correlated-diagonal theorem using an exact physical "
            "split-cubic endpoint row/column factor and explicit rational "
            "upper bounds for every remaining distinctness mask"
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
        "q64 masked cubic-endpoint repair: "
        f"repaired={result.repaired_entries},"
        f"septimic={result.cubic_septimic_entries},"
        f"recovered={result.recovered_cubic_quintic_entries},"
        f"coefficient={result.coefficient:.12g},"
        f"remaining={result.remaining_quarantined_entries}"
    )


if __name__ == "__main__":
    main()
