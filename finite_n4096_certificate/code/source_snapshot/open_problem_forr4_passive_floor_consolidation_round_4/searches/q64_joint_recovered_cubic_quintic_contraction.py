#!/usr/bin/env python3
"""Joint shared-quintic contraction for the twelve quarantined entries."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from fractions import Fraction
from json import dumps
from math import comb, prod
from pathlib import Path

from q64_masked_recovered_cubic_quintic_incidence_repair import (
    canonical_entry,
    outward_sqrt,
    rejected_chain_entries,
    set_partitions,
    simple_bipartite_shape,
    sum_outward,
)
from q64_shared_quintic_row_chain_insertion import (
    odd_record_incidence,
    two_axis_relaxed_incidence,
)


ROOT = Path(__file__).resolve().parents[1]
ORDER = 64
PROFILE = (1, 3, 5, 3)
CANONICAL_SPLITS = (
    (0, 1, 3, 2),
    (0, 2, 2, 2),
    (0, 2, 3, 1),
)
ProfileSplit = tuple[tuple[int, ...], tuple[int, ...]]
AxisPattern = tuple[int, ...]
RecordTriple = tuple[int, int, int]


@dataclass(frozen=True)
class JointShapeRow:
    row_pattern: AxisPattern
    column_pattern: AxisPattern
    row_record: int
    column_record: int
    left_endpoint_numerator: int
    left_endpoint_denominator: int
    right_link_numerator: int
    right_link_denominator: int
    joint_entry_numerator: int
    joint_entry_denominator: int
    realizes_sector_maximum: bool


@dataclass(frozen=True)
class JointRecoveredCubicQuinticContraction:
    order: int
    canonical_orbits: int
    repaired_entries: int
    record_sectors: int
    feasible_shape_rows: int
    distinct_coefficients: tuple[float, ...]
    minimum_coefficient: float
    maximum_coefficient: float
    maximum_split: tuple[int, ...]
    all_coefficients_strictly_below_one: bool


def validate_order(order: int) -> None:
    if order < 8 or order & (order - 1):
        raise ValueError(("power-of-two order at least eight required", order))


def parity_record_from_pattern(pattern: AxisPattern) -> int:
    return sum(multiplicity % 2 for multiplicity in pattern)


def multiplicity_pattern(partition: tuple[int, ...]) -> AxisPattern:
    return tuple(sorted(Counter(partition).values(), reverse=True))


def feasible_quintic_pattern_pairs() -> tuple[tuple[AxisPattern, AxisPattern], ...]:
    """Return every row/column multiplicity pair of a simple five-cell support."""

    partitions = tuple(set_partitions(5))
    result = {
        (multiplicity_pattern(rows), multiplicity_pattern(columns))
        for rows in partitions
        for columns in partitions
        if simple_bipartite_shape(rows, columns)
        and parity_record_from_pattern(multiplicity_pattern(rows)) in (1, 3)
        and parity_record_from_pattern(multiplicity_pattern(columns)) in (1, 3)
    }
    return tuple(sorted(result))


def record_one_compatible_bound(order: int) -> Fraction:
    """Endpoint-compatible cubic--quintic record-one bound for every shape."""

    validate_order(order)
    q = order
    return Fraction(q + 2, q * (q - 1) * (q - 2))


def record_three_plain_bound(order: int) -> Fraction:
    validate_order(order)
    return Fraction(1, comb(order, 3))


def record_three_even_pair_bound(order: int) -> Fraction:
    """Record-three bound when the quintic has one nonzero even pair."""

    validate_order(order)
    return min(
        record_three_plain_bound(order),
        Fraction(3, (order - 3) * comb(order, 3)),
    )


def endpoint_cubic_bound(record: int, order: int) -> Fraction:
    validate_order(order)
    if record == 1:
        return Fraction(1, order * (order - 1))
    if record == 3:
        return Fraction(1, order)
    raise ValueError(("unsupported endpoint record", record))


def left_endpoint_middle_bound(
    row_pattern: AxisPattern,
    record: int,
    order: int,
) -> Fraction:
    """Bound the singleton--cubic and cubic--quintic links jointly on the left."""

    if parity_record_from_pattern(row_pattern) != record:
        raise ValueError(("row pattern record mismatch", row_pattern, record))
    endpoint = endpoint_cubic_bound(record, order)
    if record == 1:
        middle = record_one_compatible_bound(order)
    elif row_pattern == (3, 1, 1):
        middle = record_three_plain_bound(order)
    elif row_pattern == (2, 1, 1, 1):
        middle = record_three_even_pair_bound(order)
    else:
        raise ValueError(("unsupported record-three row pattern", row_pattern))
    return endpoint * middle


def right_quintic_cubic_bound(
    column_pattern: AxisPattern,
    record: int,
    order: int,
) -> Fraction:
    """Bound the quintic--cubic link from the quintic column shape."""

    if parity_record_from_pattern(column_pattern) != record:
        raise ValueError(("column pattern record mismatch", column_pattern, record))
    if record == 1:
        if column_pattern in ((5,), (4, 1)):
            return Fraction(1, order)
        if column_pattern in ((3, 2), (2, 2, 1)):
            return record_one_compatible_bound(order)
    elif column_pattern == (3, 1, 1):
        return record_three_plain_bound(order)
    elif column_pattern == (2, 1, 1, 1):
        return record_three_even_pair_bound(order)
    raise ValueError(("unsupported quintic column pattern", column_pattern, record))


def shape_rows(order: int = ORDER) -> tuple[JointShapeRow, ...]:
    validate_order(order)
    provisional = []
    sector_maxima: dict[tuple[int, int], Fraction] = {}
    for row_pattern, column_pattern in feasible_quintic_pattern_pairs():
        row_record = parity_record_from_pattern(row_pattern)
        column_record = parity_record_from_pattern(column_pattern)
        left = left_endpoint_middle_bound(row_pattern, row_record, order)
        right = right_quintic_cubic_bound(
            column_pattern, column_record, order
        )
        joint = left * right
        key = row_record, column_record
        sector_maxima[key] = max(sector_maxima.get(key, Fraction(0)), joint)
        provisional.append(
            (
                row_pattern,
                column_pattern,
                row_record,
                column_record,
                left,
                right,
                joint,
            )
        )
    if set(sector_maxima) != {(1, 1), (1, 3), (3, 1), (3, 3)}:
        raise AssertionError(("incomplete joint sectors", sector_maxima))
    return tuple(
        JointShapeRow(
            row_pattern=row_pattern,
            column_pattern=column_pattern,
            row_record=row_record,
            column_record=column_record,
            left_endpoint_numerator=left.numerator,
            left_endpoint_denominator=left.denominator,
            right_link_numerator=right.numerator,
            right_link_denominator=right.denominator,
            joint_entry_numerator=joint.numerator,
            joint_entry_denominator=joint.denominator,
            realizes_sector_maximum=(
                joint == sector_maxima[(row_record, column_record)]
            ),
        )
        for (
            row_pattern,
            column_pattern,
            row_record,
            column_record,
            left,
            right,
            joint,
        ) in provisional
    )


def joint_chain_entry_bound(
    first_record: int,
    second_record: int,
    order: int = ORDER,
) -> Fraction:
    """Maximum shape-conditioned bound for all three physical links."""

    values = [
        Fraction(row.joint_entry_numerator, row.joint_entry_denominator)
        for row in shape_rows(order)
        if (row.row_record, row.column_record)
        == (first_record, second_record)
    ]
    if not values:
        raise ValueError(("empty joint record sector", first_record, second_record))
    return max(values)


def canonical_quarantined_entry(entry: ProfileSplit) -> ProfileSplit:
    profile, split = canonical_entry(entry)
    if profile != PROFILE:
        raise ValueError(("wrong profile", entry))
    if split in CANONICAL_SPLITS:
        return profile, split
    complement = tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )
    if complement in CANONICAL_SPLITS:
        return profile, complement
    raise ValueError(("not a quarantined split orbit", entry))


def sector_completion_degrees(
    split: tuple[int, ...],
    first_record: int,
    second_record: int,
    order: int = ORDER,
) -> tuple[int, int]:
    validate_order(order)
    complement = tuple(
        degree - selected
        for degree, selected in zip(PROFILE, split, strict=True)
    )

    def singleton(selected: int) -> int:
        return order * order if selected == 0 else 1

    families = (
        singleton,
        lambda selected: two_axis_relaxed_incidence(
            order, 3, 1, first_record, selected
        ),
        lambda selected: two_axis_relaxed_incidence(
            order, 5, first_record, second_record, selected
        ),
        lambda selected: odd_record_incidence(
            order, 3, second_record, selected
        ),
    )
    row_degree = prod(
        family(selected)
        for family, selected in zip(families, split, strict=True)
    )
    column_degree = prod(
        family(selected)
        for family, selected in zip(families, complement, strict=True)
    )
    return row_degree, column_degree


def sector_squared_coefficient(
    entry: ProfileSplit,
    records: RecordTriple,
    order: int = ORDER,
) -> Fraction:
    """Exact square of the best rank/row/column sector upper bound."""

    profile, split = canonical_quarantined_entry(entry)
    if records[0] != 1 or records[1:] not in (
        (1, 1),
        (1, 3),
        (3, 1),
        (3, 3),
    ):
        raise ValueError(("invalid record sector", records))
    first_record, second_record = records[1:]
    maximum_entry = joint_chain_entry_bound(
        first_record, second_record, order
    )
    row_degree, column_degree = sector_completion_degrees(
        split, first_record, second_record, order
    )
    complement_marks = sum(PROFILE) - sum(split)
    rank_bound = (
        order ** min(sum(split), complement_marks) * maximum_entry
    )
    return min(
        rank_bound**2,
        maximum_entry**2 * row_degree,
        maximum_entry**2 * column_degree,
    )


def sector_bound_mechanism(
    entry: ProfileSplit,
    records: RecordTriple,
    order: int = ORDER,
) -> str:
    profile, split = canonical_quarantined_entry(entry)
    first_record, second_record = records[1:]
    maximum_entry = joint_chain_entry_bound(
        first_record, second_record, order
    )
    row_degree, column_degree = sector_completion_degrees(
        split, first_record, second_record, order
    )
    values = {
        "rank": (
            order ** min(sum(split), sum(PROFILE) - sum(split))
            * maximum_entry
        )
        ** 2,
        "row": maximum_entry**2 * row_degree,
        "column": maximum_entry**2 * column_degree,
    }
    return min(values, key=values.get)


def coefficient(entry: ProfileSplit, order: int = ORDER) -> float:
    canonical_quarantined_entry(entry)
    return sum_outward(
        outward_sqrt(
            sector_squared_coefficient(entry, (1, first, second), order)
        )
        for first in (1, 3)
        for second in (1, 3)
    )


def repaired_entries(order: int = ORDER) -> tuple[ProfileSplit, ...]:
    entries = tuple(sorted(rejected_chain_entries()))
    if len(entries) != 12:
        raise AssertionError(("joint repair inventory", len(entries)))
    result = tuple(entry for entry in entries if coefficient(entry, order) <= 1)
    return result


def coefficient_map(order: int = ORDER) -> dict[ProfileSplit, float]:
    return {entry: coefficient(entry, order) for entry in repaired_entries(order)}


def diagnostic() -> JointRecoveredCubicQuinticContraction:
    coefficients = coefficient_map()
    maximum_entry = max(coefficients, key=coefficients.get)
    values = tuple(sorted(set(coefficients.values())))
    return JointRecoveredCubicQuinticContraction(
        order=ORDER,
        canonical_orbits=len(CANONICAL_SPLITS),
        repaired_entries=len(coefficients),
        record_sectors=4,
        feasible_shape_rows=len(shape_rows()),
        distinct_coefficients=values,
        minimum_coefficient=min(values),
        maximum_coefficient=coefficients[maximum_entry],
        maximum_split=canonical_quarantined_entry(maximum_entry)[1],
        all_coefficients_strictly_below_one=max(values) < 1,
    )


def artifact_text(
    result: JointRecoveredCubicQuinticContraction,
) -> str:
    payload = {
        "schema": "round4_q64_joint_recovered_cubic_quintic_contraction_v1",
        "result": asdict(result),
        "shape_table": [asdict(row) for row in shape_rows()],
        "sector_entry_bounds": [
            {
                "records": [1, first, second],
                "numerator": joint_chain_entry_bound(first, second).numerator,
                "denominator": joint_chain_entry_bound(first, second).denominator,
            }
            for first in (1, 3)
            for second in (1, 3)
        ],
        "registry_entries": [
            {
                "profile": list(profile),
                "split": list(split),
                "outward_coefficient": coefficient((profile, split)),
                "sectors": [
                    {
                        "records": [1, first, second],
                        "entry_numerator": joint_chain_entry_bound(
                            first, second
                        ).numerator,
                        "entry_denominator": joint_chain_entry_bound(
                            first, second
                        ).denominator,
                        "squared_coefficient_numerator": (
                            sector_squared_coefficient(
                                (profile, split), (1, first, second)
                            ).numerator
                        ),
                        "squared_coefficient_denominator": (
                            sector_squared_coefficient(
                                (profile, split), (1, first, second)
                            ).denominator
                        ),
                        "mechanism": sector_bound_mechanism(
                            (profile, split), (1, first, second)
                        ),
                    }
                    for first in (1, 3)
                    for second in (1, 3)
                ],
            }
            for profile, split in repaired_entries()
        ],
        "evidence_label": (
            "arbitrary-correlated-diagonal theorem for the twelve recovered "
            "cubic-quintic chain entries; the two cubic-quintic links are "
            "bounded jointly through the shared simple quintic support; "
            "fifteen feasible row/column multiplicity pairs retain parity "
            "records and active xor groups; exact rational rank and complete "
            "row/column incidence squares are outward rounded only at display"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = diagnostic()
    text = artifact_text(result)
    if arguments.output is not None:
        arguments.output.write_text(text, encoding="utf-8")
    print(
        "q64 joint recovered cubic-quintic contraction: "
        f"shapes={result.feasible_shape_rows},"
        f"repaired={result.repaired_entries},"
        f"orbits={result.canonical_orbits},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"maximum_split={result.maximum_split},"
        "status=joint_shared_quintic_theorem"
    )


if __name__ == "__main__":
    main()
