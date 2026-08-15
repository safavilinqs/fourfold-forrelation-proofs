#!/usr/bin/env python3
"""Exact generic record-chain screen for the final eighty q64 entries."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import comb, prod

from q64_degree_ten_completion_row_insertion import orbit
from q64_masked_recovered_cubic_quintic_incidence_repair import (
    block_incidence_table,
    outward_sqrt,
    sum_outward,
)
from q64_same_side_whole_link_insertion import remaining_entries


ORDER = 64
ProfileSplit = tuple[tuple[int, ...], tuple[int, ...]]
RecordTriple = tuple[int, int, int]


@dataclass(frozen=True)
class UniversalResidualScreen:
    order: int
    entries: int
    orbits: int
    profiles: int
    record_sectors: int
    passing_entries: int
    passing_orbits: int
    minimum_coefficient: float
    maximum_coefficient: float
    maximum_entry: ProfileSplit


def link_records(left_degree: int, right_degree: int) -> tuple[int, ...]:
    return tuple(range(1, min(left_degree, right_degree) + 1, 2))


def record_sectors(profile: tuple[int, ...]) -> tuple[RecordTriple, ...]:
    if len(profile) != 4:
        raise ValueError(("four-block profile required", profile))
    return tuple(
        product(
            link_records(profile[0], profile[1]),
            link_records(profile[1], profile[2]),
            link_records(profile[2], profile[3]),
        )
    )


def universal_link_bound(order: int, record: int) -> Fraction:
    return Fraction(1, comb(order, record))


def sector_entry_bound(
    profile: tuple[int, ...],
    records: RecordTriple,
    order: int = ORDER,
) -> Fraction:
    if records not in record_sectors(profile):
        raise ValueError(("invalid record sector", profile, records))
    return prod(universal_link_bound(order, record) for record in records)


def block_incidence(
    order: int,
    degree: int,
    selected: int,
    left_record: int | None,
    right_record: int | None,
) -> int:
    if left_record is None and right_record is None:
        raise ValueError("at least one adjacent link record is required")
    return block_incidence_table(
        degree,
        left_record,
        right_record,
        order,
    )[selected]


def sector_completion_degrees(
    entry: ProfileSplit,
    records: RecordTriple,
    order: int = ORDER,
) -> tuple[int, int]:
    profile, split = entry
    if records not in record_sectors(profile):
        raise ValueError(("invalid record sector", entry, records))
    complement = tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )
    adjacent = (
        (None, records[0]),
        (records[0], records[1]),
        (records[1], records[2]),
        (records[2], None),
    )

    def total(selected_sizes: tuple[int, ...]) -> int:
        return prod(
            block_incidence(
                order,
                degree,
                selected,
                left_record,
                right_record,
            )
            for degree, selected, (left_record, right_record) in zip(
                profile,
                selected_sizes,
                adjacent,
                strict=True,
            )
        )

    return total(split), total(complement)


def sector_squared_coefficient(
    entry: ProfileSplit,
    records: RecordTriple,
    order: int = ORDER,
) -> Fraction:
    profile, split = entry
    complement_size = sum(profile) - sum(split)
    entry_bound = sector_entry_bound(profile, records, order)
    row_degree, column_degree = sector_completion_degrees(
        entry, records, order
    )
    rank = order ** min(sum(split), complement_size) * entry_bound
    return min(
        rank**2,
        entry_bound**2 * row_degree,
        entry_bound**2 * column_degree,
    )


def coefficient(entry: ProfileSplit, order: int = ORDER) -> float:
    profile, _ = entry
    return sum_outward(
        outward_sqrt(sector_squared_coefficient(entry, records, order))
        for records in record_sectors(profile)
    )


def coefficient_map(order: int = ORDER) -> dict[ProfileSplit, float]:
    return {entry: coefficient(entry, order) for entry in remaining_entries()}


def diagnostic() -> UniversalResidualScreen:
    coefficients = coefficient_map()
    maximum_entry = max(coefficients, key=coefficients.get)
    passing = {entry for entry, value in coefficients.items() if value <= 1}
    orbit_keys = {min(orbit(entry)) for entry in coefficients}
    passing_orbits = {
        key
        for key in orbit_keys
        if all(entry in passing for entry in orbit(key))
    }
    return UniversalResidualScreen(
        order=ORDER,
        entries=len(coefficients),
        orbits=len(orbit_keys),
        profiles=len(Counter(profile for profile, _ in coefficients)),
        record_sectors=sum(
            len(record_sectors(profile))
            for profile in {profile for profile, _ in coefficients}
        ),
        passing_entries=len(passing),
        passing_orbits=len(passing_orbits),
        minimum_coefficient=min(coefficients.values()),
        maximum_coefficient=coefficients[maximum_entry],
        maximum_entry=maximum_entry,
    )


def main() -> None:
    result = diagnostic()
    print(
        "q64 final residual universal-chain screen: "
        f"entries={result.entries},orbits={result.orbits},"
        f"profiles={result.profiles},sectors={result.record_sectors},"
        f"passing={result.passing_entries}/{result.passing_orbits},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"maximum_entry={result.maximum_entry},"
        "status=generic_actual_mask_screen"
    )


if __name__ == "__main__":
    main()
