#!/usr/bin/env python3
"""Repair all recovered cubic--quintic entries by masked endpoint chains."""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict, dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from json import dumps
from math import comb, factorial, inf, nextafter, prod, sqrt
from pathlib import Path

from q64_masked_cubic_septimic_chain_repair import (
    coefficient_one_dependent_entries,
    cubic_endpoint_slice_energies,
    cubic_fixed_one_record_energies,
    previous_repaired_entries as pre_cubic_septimic_repaired_entries,
    repaired_entries as cubic_septimic_repaired_entries,
)
from q64_masked_double_quintic_record_repair import endpoint_record_energies
from q64_masked_four_cubic_incidence_repair import (
    record_one_link_bound,
    record_three_link_bound,
)
from q64_noncubic_recovered_universal_insertion import recovered_universal_entries
from q64_shared_quintic_row_chain_insertion import (
    odd_record_incidence,
    two_axis_relaxed_incidence,
)


ROOT = Path(__file__).resolve().parents[1]
ORDER = 64
PASSING_PROFILE = (1, 5, 3, 3)
ProfileSplit = tuple[tuple[int, ...], tuple[int, ...]]
Partition = tuple[int, ...]
RecordTriple = tuple[int, int, int]


@dataclass(frozen=True)
class MaskedRecoveredCubicQuinticIncidenceRepair:
    order: int
    candidate_entries: int
    repaired_entries: int
    repaired_orbits: int
    endpoint_refined_entries: int
    chain_refined_entries: int
    rejected_chain_entries: int
    record_sectors: int
    distinct_coefficients: tuple[float, ...]
    minimum_coefficient: float
    maximum_coefficient: float
    maximum_split: tuple[int, ...]
    remaining_quarantined_entries: int


def canonical(values: Partition) -> Partition:
    labels: dict[int, int] = {}
    return tuple(labels.setdefault(value, len(labels)) for value in values)


def set_partitions(size: int):
    """Yield restricted-growth encodings of all set partitions."""

    if size == 0:
        yield ()
        return
    values = [0] * size

    def visit(index: int, maximum: int):
        if index == size:
            yield tuple(values)
            return
        for value in range(maximum + 2):
            values[index] = value
            yield from visit(index + 1, max(maximum, value))

    values[0] = 0
    yield from visit(1, 0)


def parity_record(partition: Partition) -> int:
    counts: dict[int, int] = defaultdict(int)
    for value in partition:
        counts[value] += 1
    return sum(count % 2 for count in counts.values())


def simple_bipartite_shape(rows: Partition, columns: Partition) -> bool:
    return len(set(zip(rows, columns, strict=True))) == len(rows)


def falling(total: int, selected: int) -> int:
    return prod(range(total - selected + 1, total + 1))


@lru_cache(maxsize=None)
def block_incidence_table(
    degree: int,
    left_record: int | None,
    right_record: int | None,
    order: int = ORDER,
) -> tuple[int, ...]:
    """Exact physical completion incidences for one support block.

    ``left_record`` constrains the row parity record seen by the incoming
    plant link, and ``right_record`` constrains the column parity record seen
    by the outgoing link.  ``None`` leaves the exterior endpoint axis free.
    """

    partitions = tuple(set_partitions(degree))
    accumulators = [defaultdict(int) for _ in range(degree + 1)]
    for rows in partitions:
        if left_record is not None and parity_record(rows) != left_record:
            continue
        total_rows = len(set(rows))
        for columns in partitions:
            if right_record is not None and parity_record(columns) != right_record:
                continue
            if not simple_bipartite_shape(rows, columns):
                continue
            total_columns = len(set(columns))
            for selected in range(degree + 1):
                partial_rows = canonical(rows[:selected])
                partial_columns = canonical(columns[:selected])
                old_rows = len(set(rows[:selected]))
                old_columns = len(set(columns[:selected]))
                accumulators[selected][(partial_rows, partial_columns)] += (
                    falling(order - old_rows, total_rows - old_rows)
                    * falling(
                        order - old_columns,
                        total_columns - old_columns,
                    )
                )
    result = []
    for selected, values in enumerate(accumulators):
        divisor = factorial(degree - selected)
        if not values or any(value % divisor for value in values.values()):
            raise AssertionError(
                ("invalid block incidence", degree, left_record, right_record)
            )
        result.append(max(value // divisor for value in values.values()))
    return tuple(result)


def link_records(left_degree: int, right_degree: int) -> tuple[int, ...]:
    return tuple(range(1, min(left_degree, right_degree) + 1, 2))


def link_bound(
    left_degree: int,
    right_degree: int,
    record: int,
    order: int = ORDER,
) -> Fraction:
    """Uniform one-link moment bound in a fixed common-record sector."""

    if record not in link_records(left_degree, right_degree):
        raise ValueError(("incompatible link record", left_degree, right_degree, record))
    if left_degree == right_degree == 3:
        # The smaller record-one formula is not universal: a vertical
        # cubic followed by a horizontal cubic has moment 1/q.  It may only
        # be used after an adjacent record-one constraint excludes that
        # geometry.  This context-free helper therefore keeps the generic
        # signed-permutation bound.
        if record == 1:
            return Fraction(1, order)
        return record_three_link_bound(order)
    return Fraction(1, comb(order, record))


def compatible_cubic_cubic_record_one_bound(
    order: int = ORDER,
) -> Fraction:
    """Record-one cubic--cubic bound with a record-one exterior axis."""

    return record_one_link_bound(order)


def endpoint_remaining_link_bound(
    profile: tuple[int, ...],
    records: RecordTriple,
    link_index: int,
    order: int = ORDER,
) -> Fraction:
    """Safe remaining-link bound for the endpoint-row mechanism.

    The improved cubic--cubic record-one estimate requires the left cubic's
    incoming record also to be one.  In every other record-one cubic--cubic
    sector the generic coefficient 1/q is used.
    """

    left_degree = profile[link_index]
    right_degree = profile[link_index + 1]
    record = records[link_index]
    if left_degree == right_degree == 3 and record == 1:
        incoming_record = records[link_index - 1]
        if incoming_record == 1:
            return compatible_cubic_cubic_record_one_bound(order)
        return Fraction(1, order)
    return link_bound(left_degree, right_degree, record, order)


def record_sectors(profile: tuple[int, ...]) -> tuple[RecordTriple, ...]:
    return tuple(
        product(
            *(
                link_records(left, right)
                for left, right in zip(profile, profile[1:])
            )
        )
    )


def sector_incidences(
    profile: tuple[int, ...],
    records: RecordTriple,
    order: int = ORDER,
) -> tuple[tuple[int, ...], ...]:
    result = []
    for index, degree in enumerate(profile):
        left_record = records[index - 1] if index else None
        right_record = records[index] if index < len(records) else None
        result.append(
            block_incidence_table(
                degree,
                left_record,
                right_record,
                order,
            )
        )
    return tuple(result)


def sector_maximum_entry(
    profile: tuple[int, ...],
    records: RecordTriple,
    order: int = ORDER,
) -> Fraction:
    return prod(
        link_bound(left, right, record, order)
        for left, right, record in zip(
            profile,
            profile[1:],
            records,
        )
    )


def sector_squared_coefficient(
    profile: tuple[int, ...],
    split: tuple[int, ...],
    records: RecordTriple,
    order: int = ORDER,
) -> Fraction:
    families = sector_incidences(profile, records, order)
    row_degree = prod(
        family[selected]
        for family, selected in zip(families, split, strict=True)
    )
    column_degree = prod(
        family[degree - selected]
        for family, degree, selected in zip(
            families,
            profile,
            split,
            strict=True,
        )
    )
    maximum_entry = sector_maximum_entry(profile, records, order)
    return maximum_entry**2 * min(row_degree, column_degree)


def outward_sqrt(exact: Fraction) -> float:
    value = sqrt(float(exact))
    while Fraction.from_float(value) ** 2 < exact:
        value = nextafter(value, inf)
    return value


def canonical_entry(entry: ProfileSplit) -> ProfileSplit:
    profile, split = entry
    reversed_profile = tuple(reversed(profile))
    if reversed_profile < profile:
        return reversed_profile, tuple(reversed(split))
    return profile, split


def incidence_coefficient(entry: ProfileSplit, order: int = ORDER) -> float:
    profile, split = canonical_entry(entry)
    result = 0.0
    for records in record_sectors(profile):
        term = outward_sqrt(
            sector_squared_coefficient(profile, split, records, order)
        )
        result = nextafter(result + term, inf)
    return result


def endpoint_energy(
    degree: int,
    selected: int,
    singleton_selected: int,
    outgoing_record: int,
    order: int = ORDER,
) -> Fraction:
    """Complete singleton--cubic/quintic row energy in one next-link record."""

    variable_singleton = singleton_selected == 0
    if singleton_selected not in (0, 1):
        raise ValueError(("invalid singleton split", singleton_selected))
    if degree == 5:
        energies = endpoint_record_energies(
            selected, variable_singleton, order
        )
        if energies is None:
            raise ValueError(("uncalibrated quintic endpoint", selected))
        return energies[outgoing_record]
    if degree == 3:
        multiplier = order**2 if variable_singleton else 1
        if selected == 1:
            return (
                multiplier
                * cubic_fixed_one_record_energies(order)[outgoing_record]
            )
        return multiplier * cubic_endpoint_slice_energies(order)[selected]
    raise ValueError(("unsupported endpoint degree", degree))


def endpoint_row_sector_squared(
    profile: tuple[int, ...],
    split: tuple[int, ...],
    records: RecordTriple,
    order: int = ORDER,
) -> Fraction:
    """Exact-energy endpoint factor followed by two physical block rows."""

    if profile[0] != 1:
        raise ValueError(("canonical endpoint singleton required", profile))
    energy = endpoint_energy(
        profile[1], split[1], split[0], records[1], order
    )
    families = sector_incidences(profile, records, order)
    remaining_incidence = prod(
        families[index][split[index]] for index in (2, 3)
    )
    remaining_moment = prod(
        endpoint_remaining_link_bound(profile, records, index, order)
        for index in (1, 2)
    )
    return energy * remaining_incidence * remaining_moment**2


def endpoint_sector_squared_coefficient(
    entry: ProfileSplit,
    records: RecordTriple,
    order: int = ORDER,
) -> Fraction:
    profile, split = canonical_entry(entry)
    complement = tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )
    return min(
        endpoint_row_sector_squared(profile, split, records, order),
        endpoint_row_sector_squared(profile, complement, records, order),
    )


def endpoint_coefficient(entry: ProfileSplit, order: int = ORDER) -> float:
    profile, _ = canonical_entry(entry)
    return sum_outward(
        outward_sqrt(
            endpoint_sector_squared_coefficient(entry, records, order)
        )
        for records in record_sectors(profile)
    )


def cubic_quintic_record_bound(
    record: int, order: int = ORDER
) -> Fraction:
    """Sharp analytic cubic--quintic fixed-record moment bound."""

    q = order
    if record == 1:
        return Fraction(q + 2, q * (q - 1) * (q - 2))
    if record == 3:
        return Fraction(3, (q - 3) * comb(q, 3))
    raise ValueError(("unsupported cubic-quintic record", record))


def endpoint_cubic_record_bound(
    outgoing_record: int, order: int = ORDER
) -> Fraction:
    """Joint singleton--cubic endpoint bound in the outgoing record."""

    if outgoing_record == 1:
        return Fraction(1, order * (order - 1))
    if outgoing_record == 3:
        return Fraction(1, order)
    raise ValueError(("unsupported endpoint-cubic record", outgoing_record))


def chain_sector_squared_coefficient(
    entry: ProfileSplit,
    records: RecordTriple,
    order: int = ORDER,
) -> Fraction:
    """Withdrawn pre-audit chain expression; this is not a valid upper bound."""

    profile, split = canonical_entry(entry)
    if profile != (1, 3, 5, 3):
        raise ValueError(("not an endpoint-cubic/quintic chain", entry))
    _, first_record, second_record = records
    q = order

    def singleton(selected: int) -> int:
        return q * q if selected == 0 else 1

    families = (
        singleton,
        lambda selected: two_axis_relaxed_incidence(
            q, 3, 1, first_record, selected
        ),
        lambda selected: two_axis_relaxed_incidence(
            q, 5, first_record, second_record, selected
        ),
        lambda selected: odd_record_incidence(
            q, 3, second_record, selected
        ),
    )
    complement = tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )
    row_degree = prod(
        family(selected)
        for family, selected in zip(families, split, strict=True)
    )
    column_degree = prod(
        family(selected)
        for family, selected in zip(families, complement, strict=True)
    )
    maximum_entry = (
        endpoint_cubic_record_bound(first_record, q)
        * cubic_quintic_record_bound(first_record, q)
        * cubic_quintic_record_bound(second_record, q)
    )
    rank_bound = (
        q ** min(sum(split), sum(complement)) * maximum_entry
    )
    return min(
        rank_bound**2,
        maximum_entry**2 * row_degree,
        maximum_entry**2 * column_degree,
    )


def chain_coefficient(entry: ProfileSplit, order: int = ORDER) -> float:
    """Return the withdrawn chain expression for historical diagnostics only."""

    return sum_outward(
        outward_sqrt(
            chain_sector_squared_coefficient(entry, records, order)
        )
        for records in record_sectors(canonical_entry(entry)[0])
    )


def sum_outward(values) -> float:
    result = 0.0
    for value in values:
        result = nextafter(result + value, inf)
    return result


def coefficient_with_mechanism(
    entry: ProfileSplit, order: int = ORDER
) -> tuple[float, str]:
    profile, _ = canonical_entry(entry)
    if profile == (1, 3, 5, 3):
        return inf, "rejected_four_sector_physical_chain"
    return endpoint_coefficient(entry, order), "record_resolved_endpoint_row"


def coefficient(entry: ProfileSplit, order: int = ORDER) -> float:
    return coefficient_with_mechanism(entry, order)[0]


def previous_repaired_entries() -> frozenset[ProfileSplit]:
    return frozenset(pre_cubic_septimic_repaired_entries()) | frozenset(
        cubic_septimic_repaired_entries()
    )


def candidate_entries() -> tuple[ProfileSplit, ...]:
    previous = previous_repaired_entries()
    return tuple(
        sorted(
            entry
            for entry in recovered_universal_entries()
            if entry not in previous
        )
    )


def repaired_entries() -> tuple[ProfileSplit, ...]:
    return tuple(
        entry
        for entry in candidate_entries()
        if coefficient_with_mechanism(entry)[1]
        == "record_resolved_endpoint_row"
        and coefficient(entry) <= 1
    )


def rejected_chain_entries() -> tuple[ProfileSplit, ...]:
    return tuple(
        entry
        for entry in candidate_entries()
        if canonical_entry(entry)[0] == (1, 3, 5, 3)
    )


def coefficient_map() -> dict[ProfileSplit, float]:
    return {entry: coefficient(entry) for entry in repaired_entries()}


def diagnostic() -> MaskedRecoveredCubicQuinticIncidenceRepair:
    candidates = candidate_entries()
    repaired = repaired_entries()
    rejected = rejected_chain_entries()
    if len(candidates) != 40 or len(repaired) != 28 or len(rejected) != 12:
        raise AssertionError(("recovered incidence inventory", len(candidates), len(repaired)))
    values = {entry: coefficient(entry) for entry in repaired}
    maximum_entry = max(values, key=values.get)
    distinct = tuple(sorted(set(values.values())))
    return MaskedRecoveredCubicQuinticIncidenceRepair(
        order=ORDER,
        candidate_entries=len(candidates),
        repaired_entries=len(repaired),
        repaired_orbits=7,
        endpoint_refined_entries=sum(
            coefficient_with_mechanism(entry)[1]
            == "record_resolved_endpoint_row"
            for entry in repaired
        ),
        chain_refined_entries=0,
        rejected_chain_entries=len(rejected),
        record_sectors=len(record_sectors(PASSING_PROFILE)),
        distinct_coefficients=distinct,
        minimum_coefficient=min(values.values()),
        maximum_coefficient=values[maximum_entry],
        maximum_split=maximum_entry[1],
        remaining_quarantined_entries=(
            len(coefficient_one_dependent_entries())
            - len(previous_repaired_entries())
            - len(repaired)
        ),
    )


def artifact_text(
    result: MaskedRecoveredCubicQuinticIncidenceRepair,
) -> str:
    payload = {
        "schema": "round4_q64_masked_recovered_cubic_quintic_incidence_repair_v3",
        "result": asdict(result),
        "passing_profile_incidence_tables": {
            f"degree_{degree}_left_{left}_right_{right}": list(
                block_incidence_table(degree, left, right)
            )
            for degree, left, right in (
                (1, None, 1),
                (5, 1, 1),
                (5, 1, 3),
                (3, 1, 1),
                (3, 1, 3),
                (3, 3, 1),
                (3, 3, 3),
                (3, 1, None),
                (3, 3, None),
            )
        },
        "repaired_registry_entries": [
            {
                "profile": list(profile),
                "split": list(split),
                "outward_coefficient": coefficient((profile, split)),
                "mechanism": coefficient_with_mechanism((profile, split))[1],
                "sector_squared_coefficients": [
                    {
                        "records": list(records),
                        "numerator": endpoint_sector_squared_coefficient(
                            (profile, split), records
                        ).numerator,
                        "denominator": endpoint_sector_squared_coefficient(
                            (profile, split), records
                        ).denominator,
                    }
                    for records in record_sectors(
                        canonical_entry((profile, split))[0]
                    )
                ],
            }
            for profile, split in repaired_entries()
        ],
        "rejected_chain_entries": [
            {"profile": list(profile), "split": list(split)}
            for profile, split in rejected_chain_entries()
        ],
        "counterexample": {
            "order": 8,
            "records": [1, 3, 1],
            "supports": [[0], [0, 1, 2], [0, 1, 8, 16, 24], [0, 1, 2]],
            "exact_entry_numerator": 1,
            "exact_entry_denominator": 17920,
            "claimed_maximum_entry_numerator": 1,
            "claimed_maximum_entry_denominator": 25088,
        },
        "evidence_label": (
            "arbitrary-correlated-diagonal theorem for twenty-eight recovered "
            "cubic-quintic endpoint-row entries; twelve four-sector chain "
            "entries are rejected because an exact q8 physical entry exceeds "
            "the claimed product maximum; physical completion incidences retain "
            "every occurrence mask and all displayed square roots are outward rounded"
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
        "q64 masked recovered cubic-quintic incidence repair: "
        f"repaired={result.repaired_entries},"
        f"orbits={result.repaired_orbits},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"remaining={result.remaining_quarantined_entries}"
    )


if __name__ == "__main__":
    main()
