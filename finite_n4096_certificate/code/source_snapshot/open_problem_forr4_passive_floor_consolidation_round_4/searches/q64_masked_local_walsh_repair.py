#!/usr/bin/env python3
"""Repair 180 masked entries with a local q^-1 Walsh-chain gain.

The inherited coefficient-one argument is used only for completed unmasked
link kernels.  Every physical split-block mask is restored explicitly by an
inclusion--exclusion gamma_2 factor.  An internal singleton between two
higher blocks, or a same-side adjacent singleton pair, supplies an exact
factor 1/q.  Integer upper bounds keep the resulting q=64 coefficients
rational and outward-safe.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from json import dumps
from math import comb, isqrt
from pathlib import Path

from q64_masked_quintic_slice_repair import (
    coefficient_one_dependent_entries,
    repaired_entries as quintic_repaired_entries,
)


ROOT = Path(__file__).resolve().parents[1]
ORDER = 64
Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class MaskedLocalWalshRepair:
    order: int
    affected_entries: int
    previous_repaired_entries: int
    candidate_entries: int
    repaired_entries: int
    internal_singleton_entries: int
    same_side_singleton_pair_entries: int
    singleton_pair_chain_entries: int
    remaining_quarantined_entries: int
    distinct_coefficients: tuple[float, ...]
    minimum_coefficient: float
    maximum_coefficient: float


def ceil_sqrt(value: int) -> int:
    root = isqrt(value)
    return root if root * root == value else root + 1


def integer_cross_mask_factor(left: int, right: int) -> int:
    """Integer upper bound for a disjoint r-set/s-set Schur factor."""

    return sum(
        ceil_sqrt(comb(left, size) * comb(right, size))
        for size in range(min(left, right) + 1)
    )


def entry_integer_mask_factor(entry: ProfileSplit) -> int:
    profile, split = entry
    factor = 1
    for degree, selected in zip(profile, split, strict=True):
        if 0 < selected < degree:
            factor *= integer_cross_mask_factor(selected, degree - selected)
    return factor


def has_internal_singleton(profile: Profile) -> bool:
    return any(
        profile[index] == 1
        and profile[index - 1] > 1
        and profile[index + 1] > 1
        for index in (1, 2)
    )


def has_any_internal_singleton(profile: Profile) -> bool:
    return any(profile[index] == 1 for index in (1, 2))


def has_same_side_singleton_pair(profile: Profile, split: Split) -> bool:
    return any(
        profile[index : index + 2] == (1, 1)
        and split[index] == split[index + 1]
        for index in range(3)
    )


def mechanism(entry: ProfileSplit) -> str | None:
    profile, split = entry
    internal = has_internal_singleton(profile)
    pair = has_same_side_singleton_pair(profile, split)
    if internal:
        return "internal_singleton"
    if pair:
        return "same_side_singleton_pair"
    if has_any_internal_singleton(profile):
        return "singleton_pair_chain"
    return None


def candidate_entries() -> tuple[ProfileSplit, ...]:
    previous = set(quintic_repaired_entries())
    return tuple(
        sorted(
            entry
            for entry in coefficient_one_dependent_entries()
            if entry not in previous and mechanism(entry) is not None
        )
    )


def coefficient(entry: ProfileSplit) -> Fraction:
    if mechanism(entry) is None:
        raise ValueError(("entry has no local Walsh gain", entry))
    return Fraction(entry_integer_mask_factor(entry), ORDER)


def repaired_entries() -> tuple[ProfileSplit, ...]:
    return tuple(entry for entry in candidate_entries() if coefficient(entry) <= 1)


def coefficient_map() -> dict[ProfileSplit, float]:
    return {entry: float(coefficient(entry)) for entry in repaired_entries()}


def diagnostic() -> MaskedLocalWalshRepair:
    affected = coefficient_one_dependent_entries()
    previous = quintic_repaired_entries()
    candidates = candidate_entries()
    repaired = repaired_entries()
    values = tuple(sorted({coefficient(entry) for entry in repaired}))
    if len(candidates) != 180 or len(repaired) != 180:
        raise AssertionError(("local Walsh inventory", len(candidates), len(repaired)))
    return MaskedLocalWalshRepair(
        order=ORDER,
        affected_entries=len(affected),
        previous_repaired_entries=len(previous),
        candidate_entries=len(candidates),
        repaired_entries=len(repaired),
        internal_singleton_entries=sum(
            has_internal_singleton(entry[0]) for entry in repaired
        ),
        same_side_singleton_pair_entries=sum(
            mechanism(entry) == "same_side_singleton_pair" for entry in repaired
        ),
        singleton_pair_chain_entries=sum(
            mechanism(entry) == "singleton_pair_chain" for entry in repaired
        ),
        remaining_quarantined_entries=len(affected) - len(previous) - len(repaired),
        distinct_coefficients=tuple(map(float, values)),
        minimum_coefficient=float(min(values)),
        maximum_coefficient=float(max(values)),
    )


def artifact_text(result: MaskedLocalWalshRepair) -> str:
    payload = {
        "schema": "round4_q64_masked_local_walsh_repair_v1",
        "result": asdict(result),
        "repaired_registry_entries": [
            {
                "profile": list(profile),
                "split": list(split),
                "mechanism": mechanism((profile, split)),
                "integer_mask_factor": entry_integer_mask_factor((profile, split)),
                "coefficient_numerator": coefficient((profile, split)).numerator,
                "coefficient_denominator": coefficient((profile, split)).denominator,
            }
            for profile, split in repaired_entries()
        ],
        "evidence_label": (
            "arbitrary-correlated-diagonal theorem for actual masked kernels; "
            "one exact local q^-1 Walsh gain and integer-certified inclusion-"
            "exclusion factors for every physical split-block mask"
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
        "q64 masked local-Walsh repair: "
        f"candidates={result.candidate_entries},"
        f"repaired={result.repaired_entries},"
        f"internal={result.internal_singleton_entries},"
        f"same_side_pair={result.same_side_singleton_pair_entries},"
        f"remaining={result.remaining_quarantined_entries},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g}"
    )


if __name__ == "__main__":
    main()
