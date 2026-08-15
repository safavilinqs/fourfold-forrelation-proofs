#!/usr/bin/env python3
"""Regression for the retained endpoint repair and chain quarantine."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_degree_ten_completion_row_insertion import orbit  # noqa: E402
from q64_masked_four_cubic_incidence_repair import (  # noqa: E402
    endpoint_record_one_incidence,
    endpoint_record_three_incidence,
    middle_record_one_one_incidence,
    middle_record_one_three_incidence,
    middle_record_three_three_incidence,
)
from q64_masked_recovered_cubic_quintic_incidence_repair import (  # noqa: E402
    artifact_text,
    block_incidence_table,
    candidate_entries,
    canonical_entry,
    coefficient,
    coefficient_with_mechanism,
    coefficient_map,
    diagnostic,
    endpoint_coefficient,
    endpoint_energy,
    endpoint_remaining_link_bound,
    endpoint_sector_squared_coefficient,
    incidence_coefficient,
    link_bound,
    outward_sqrt,
    record_sectors,
    rejected_chain_entries,
    repaired_entries,
    sector_squared_coefficient,
)
from q64_shared_quintic_row_chain_insertion import (  # noqa: E402
    nonfavorable_adjacent_split_coefficient,
)


def support_records(order: int, support: tuple[int, ...]) -> tuple[int, int]:
    rows = Counter(coordinate // order for coordinate in support)
    columns = Counter(coordinate % order for coordinate in support)
    return (
        sum(count % 2 for count in rows.values()),
        sum(count % 2 for count in columns.values()),
    )


def family_incidence_maxima(
    order: int,
    degree: int,
    left_record: int | None,
    right_record: int | None,
) -> tuple[int, ...]:
    supports = tuple(combinations(range(order * order), degree))
    family = tuple(
        support
        for support in supports
        if (
            left_record is None
            or support_records(order, support)[0] == left_record
        )
        and (
            right_record is None
            or support_records(order, support)[1] == right_record
        )
    )
    result = []
    for selected in range(degree + 1):
        counts: Counter[tuple[int, ...]] = Counter()
        for support in family:
            counts.update(combinations(support, selected))
        result.append(max(counts.values()))
    return tuple(result)


def exact_q4_incidence_check() -> int:
    signatures = (
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
    for degree, left, right in signatures:
        direct = family_incidence_maxima(4, degree, left, right)
        formula = block_incidence_table(degree, left, right, 4)
        if direct != formula:
            raise AssertionError(
                ("q4 block incidence", degree, left, right, direct, formula)
            )
    return len(signatures)


def inherited_cubic_formula_check() -> int:
    for order in (4, 8, 64):
        expected = {
            (None, 1): endpoint_record_one_incidence(order),
            (None, 3): endpoint_record_three_incidence(order),
            (1, None): endpoint_record_one_incidence(order),
            (3, None): endpoint_record_three_incidence(order),
            (1, 1): middle_record_one_one_incidence(order),
            (1, 3): middle_record_one_three_incidence(order),
            (3, 1): middle_record_one_three_incidence(order),
            (3, 3): middle_record_three_three_incidence(order),
        }
        for records, value in expected.items():
            observed = block_incidence_table(3, *records, order)
            if observed != value:
                raise AssertionError(
                    ("inherited cubic incidence", order, records, observed, value)
                )
    return 3 * len(expected)


def main() -> None:
    candidates = candidate_entries()
    repaired = repaired_entries()
    repaired_set = set(repaired)
    rejected = rejected_chain_entries()
    if len(candidates) != 40 or len(repaired) != 28 or len(rejected) != 12:
        raise AssertionError(("recovered incidence inventory", len(candidates), len(repaired)))
    if any(not set(orbit(entry)).issubset(repaired_set) for entry in repaired):
        raise AssertionError("recovered incidence repair is not orbit closed")
    if len({frozenset(orbit(entry)) for entry in repaired}) != 7:
        raise AssertionError("recovered incidence orbit count")
    if set(coefficient_map()) != repaired_set:
        raise AssertionError("recovered incidence coefficient map")

    expected_coefficients = (
        0.016288857182016264,
        0.026491319562908268,
        0.033190827461608174,
        0.05064816804156321,
        0.11089844573555295,
        0.41391074302652664,
        0.7036151810879134,
    )
    observed_coefficients = tuple(sorted(set(coefficient_map().values())))
    if observed_coefficients != expected_coefficients:
        raise AssertionError(("recovered incidence coefficients", observed_coefficients))
    mechanisms = Counter(
        coefficient_with_mechanism(entry)[1] for entry in repaired
    )
    if mechanisms != {"record_resolved_endpoint_row": 28}:
        raise AssertionError(("recovered mechanisms", mechanisms))

    for entry in repaired:
        profile, split = entry
        canonical_profile, canonical_split = canonical_entry(entry)
        for records in record_sectors(canonical_profile):
            exact = endpoint_sector_squared_coefficient(entry, records)
            # The displayed coefficient is a sum of individually outward
            # square-rooted exact sector values.
            if Fraction.from_float(outward_sqrt(exact)) ** 2 < exact:
                raise AssertionError(("sector square root not outward", profile, split, records))
        if coefficient(entry) > 1:
            raise AssertionError(("recovered incidence coefficient above one", profile, split))
        if endpoint_coefficient(entry) > incidence_coefficient(entry):
            raise AssertionError(("endpoint refinement did not improve", entry))

    if any(coefficient(entry) != float("inf") for entry in rejected):
        raise AssertionError("rejected chain entry received a theorem coefficient")
    if any(
        coefficient_with_mechanism(entry)[1]
        != "rejected_four_sector_physical_chain"
        for entry in rejected
    ):
        raise AssertionError("rejected chain mechanism label")

    for record in (1, 3):
        if endpoint_energy(3, 1, 0, record, 4) != (
            16 * endpoint_energy(3, 1, 1, record, 4)
        ):
            raise AssertionError(("variable cubic endpoint transitivity", record))

    if link_bound(5, 3, 1) != Fraction(1, 64):
        raise AssertionError("generic record-one link bound")
    if link_bound(5, 3, 3) != Fraction(1, 41664):
        raise AssertionError("generic record-three link bound")
    if endpoint_remaining_link_bound(
        (1, 5, 3, 3), (1, 3, 1), 2
    ) != Fraction(1, 64):
        raise AssertionError("unsafe cubic record-one refinement")

    q4_families = exact_q4_incidence_check()
    cubic_formula_checks = inherited_cubic_formula_check()
    result = diagnostic()
    if (
        result.repaired_entries,
        result.repaired_orbits,
        result.record_sectors,
        result.maximum_split,
        result.remaining_quarantined_entries,
        result.endpoint_refined_entries,
        result.chain_refined_entries,
        result.rejected_chain_entries,
    ) != (28, 7, 4, (0, 1, 2, 3), 24, 28, 0, 12):
        raise AssertionError(("recovered incidence diagnostic", result))
    committed = (
        ROOT
        / "artifacts"
        / "q64_masked_recovered_cubic_quintic_incidence_repair.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale recovered incidence artifact")
    print(
        "q64 masked recovered cubic-quintic incidence regression passed: "
        f"repaired={result.repaired_entries},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"remaining={result.remaining_quarantined_entries},"
        f"endpoint={result.endpoint_refined_entries},"
        f"rejected_chain={result.rejected_chain_entries},"
        f"q4_families={q4_families},"
        f"cubic_formula_checks={cubic_formula_checks}"
    )


if __name__ == "__main__":
    main()
