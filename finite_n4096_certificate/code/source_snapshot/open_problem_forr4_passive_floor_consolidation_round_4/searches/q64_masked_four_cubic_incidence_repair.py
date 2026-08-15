#!/usr/bin/env python3
"""Repair all 38 residual four-cubic entries by masked incidences."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from itertools import product
from json import dumps
from math import comb, inf, nextafter, prod, sqrt
from pathlib import Path

from q64_masked_cubic_endpoint_repair import (
    repaired_entries as cubic_endpoint_repaired_entries,
)
from q64_masked_double_quintic_endpoint_repair import (
    repaired_entries as double_quintic_endpoint_repaired_entries,
)
from q64_masked_local_walsh_repair import (
    repaired_entries as local_walsh_repaired_entries,
)
from q64_masked_quintic_slice_repair import (
    coefficient_one_dependent_entries,
    repaired_entries as quintic_repaired_entries,
)
from q64_universal_double_cubic_insertion import double_cubic_entries
from q64_universal_multicubic_insertion import multicubic_entries


ROOT = Path(__file__).resolve().parents[1]
ORDER = 64
PROFILE = (3, 3, 3, 3)
RECORDS = (1, 3)
ProfileSplit = tuple[tuple[int, ...], tuple[int, ...]]
RecordTriple = tuple[int, int, int]


@dataclass(frozen=True)
class MaskedFourCubicIncidenceRepair:
    order: int
    repaired_entries: int
    multicubic_entries: int
    double_cubic_entries: int
    record_sectors: int
    minimum_coefficient: float
    maximum_coefficient: float
    maximum_split: tuple[int, ...]
    remaining_quarantined_entries: int


def previous_repaired_entries() -> frozenset[ProfileSplit]:
    return (
        frozenset(quintic_repaired_entries())
        | frozenset(local_walsh_repaired_entries())
        | frozenset(cubic_endpoint_repaired_entries())
        | frozenset(double_quintic_endpoint_repaired_entries())
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


def record_one_link_bound(order: int = ORDER) -> Fraction:
    q = order
    return Fraction(q + 2, q * (q - 1) * (q - 2))


def record_three_link_bound(order: int = ORDER) -> Fraction:
    return Fraction(1, comb(order, 3))


def endpoint_record_one_incidence(order: int = ORDER) -> tuple[int, ...]:
    q = order
    return (
        q * comb(q, 3) + q * q * (q - 1) * comb(q, 2),
        comb(q - 1, 2) + (q - 1) * comb(q, 2) + q * (q - 1) ** 2,
        q * q - 2,
        1,
    )


def endpoint_record_three_incidence(order: int = ORDER) -> tuple[int, ...]:
    q = order
    return (
        comb(q, 3) * q**3,
        comb(q - 1, 2) * q**2,
        q * (q - 2),
        1,
    )


def middle_record_one_one_incidence(order: int = ORDER) -> tuple[int, ...]:
    q = order
    return (q * q * (q - 1) ** 2, 3 * (q - 1) ** 2, 2 * (q - 1), 1)


def middle_record_one_three_incidence(order: int = ORDER) -> tuple[int, ...]:
    q = order
    return (
        q * comb(q, 3) * (3 * q - 2),
        (q - 1) * (q - 2) * (3 * q - 2) // 2,
        q * (q - 2),
        1,
    )


def middle_record_three_three_incidence(order: int = ORDER) -> tuple[int, ...]:
    q = order
    return (
        6 * comb(q, 3) ** 2,
        2 * comb(q - 1, 2) ** 2,
        (q - 2) ** 2,
        1,
    )


def endpoint_incidence(record: int, order: int = ORDER) -> tuple[int, ...]:
    if record == 1:
        return endpoint_record_one_incidence(order)
    if record == 3:
        return endpoint_record_three_incidence(order)
    raise ValueError(("cubic record must be one or three", record))


def middle_incidence(
    left_record: int,
    right_record: int,
    order: int = ORDER,
) -> tuple[int, ...]:
    pair = (left_record, right_record)
    if pair == (1, 1):
        return middle_record_one_one_incidence(order)
    if pair in ((1, 3), (3, 1)):
        return middle_record_one_three_incidence(order)
    if pair == (3, 3):
        return middle_record_three_three_incidence(order)
    raise ValueError(("cubic records must be one or three", pair))


def sector_maximum_entry(
    records: RecordTriple,
    order: int = ORDER,
) -> Fraction:
    bounds = {
        1: record_one_link_bound(order),
        3: record_three_link_bound(order),
    }
    return prod(bounds[record] for record in records)


def sector_incidences(
    records: RecordTriple,
    order: int = ORDER,
) -> tuple[tuple[int, ...], ...]:
    first, second, third = records
    return (
        endpoint_incidence(first, order),
        middle_incidence(first, second, order),
        middle_incidence(second, third, order),
        endpoint_incidence(third, order),
    )


def sector_squared_coefficient(
    split: tuple[int, ...],
    records: RecordTriple,
    order: int = ORDER,
) -> Fraction:
    families = sector_incidences(records, order)
    row_degree = prod(
        family[selected]
        for family, selected in zip(families, split, strict=True)
    )
    column_degree = prod(
        family[3 - selected]
        for family, selected in zip(families, split, strict=True)
    )
    maximum_entry = sector_maximum_entry(records, order)
    return maximum_entry**2 * min(row_degree, column_degree)


def outward_sqrt(exact: Fraction) -> float:
    value = sqrt(float(exact))
    while Fraction.from_float(value) ** 2 < exact:
        value = nextafter(value, inf)
    return value


def coefficient(entry: ProfileSplit, order: int = ORDER) -> float:
    profile, split = entry
    if profile != PROFILE:
        raise ValueError(("not a four-cubic entry", entry))
    result = 0.0
    for records in product(RECORDS, repeat=3):
        term = outward_sqrt(sector_squared_coefficient(split, records, order))
        result = nextafter(result + term, inf)
    return result


def coefficient_map() -> dict[ProfileSplit, float]:
    return {entry: coefficient(entry) for entry in candidate_entries()}


def repaired_entries() -> tuple[ProfileSplit, ...]:
    return tuple(entry for entry in candidate_entries() if coefficient(entry) <= 1)


def diagnostic() -> MaskedFourCubicIncidenceRepair:
    candidates = candidate_entries()
    repaired = repaired_entries()
    multicubic = frozenset(multicubic_entries())
    double_cubic = frozenset(double_cubic_entries())
    if len(candidates) != 38 or len(repaired) != 38:
        raise AssertionError(("four-cubic inventory", len(candidates), len(repaired)))
    values = {entry: coefficient(entry) for entry in repaired}
    maximum_entry = max(values, key=values.get)
    return MaskedFourCubicIncidenceRepair(
        order=ORDER,
        repaired_entries=len(repaired),
        multicubic_entries=len(set(repaired).intersection(multicubic)),
        double_cubic_entries=len(set(repaired).intersection(double_cubic)),
        record_sectors=len(tuple(product(RECORDS, repeat=3))),
        minimum_coefficient=min(values.values()),
        maximum_coefficient=values[maximum_entry],
        maximum_split=maximum_entry[1],
        remaining_quarantined_entries=(
            len(coefficient_one_dependent_entries())
            - len(previous_repaired_entries())
            - len(repaired)
        ),
    )


def artifact_text(result: MaskedFourCubicIncidenceRepair) -> str:
    payload = {
        "schema": "round4_q64_masked_four_cubic_incidence_repair_v1",
        "result": asdict(result),
        "repaired_registry_entries": [
            {
                "profile": list(profile),
                "split": list(split),
                "outward_coefficient": coefficient((profile, split)),
                "sector_squared_coefficients": [
                    {
                        "records": list(records),
                        "numerator": sector_squared_coefficient(split, records).numerator,
                        "denominator": sector_squared_coefficient(split, records).denominator,
                    }
                    for records in product(RECORDS, repeat=3)
                ],
            }
            for profile, split in repaired_entries()
        ],
        "evidence_label": (
            "arbitrary-correlated-diagonal theorem using eight exact cubic "
            "record sectors, exact physical completion-incidence bounds, "
            "and outward-rounded row/column feature coefficients"
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
        "q64 masked four-cubic incidence repair: "
        f"repaired={result.repaired_entries},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"remaining={result.remaining_quarantined_entries}"
    )


if __name__ == "__main__":
    main()
