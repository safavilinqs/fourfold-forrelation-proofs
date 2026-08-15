#!/usr/bin/env python3
"""Regression for the 38-entry physical four-cubic incidence repair."""

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
    artifact_text,
    candidate_entries,
    coefficient,
    coefficient_map,
    diagnostic,
    endpoint_record_one_incidence,
    endpoint_record_three_incidence,
    middle_record_one_one_incidence,
    middle_record_one_three_incidence,
    middle_record_three_three_incidence,
    outward_sqrt,
    record_one_link_bound,
    record_three_link_bound,
    repaired_entries,
    sector_squared_coefficient,
)
from signed_permutation_link_moment import moment  # noqa: E402


def support_records(order: int, support: tuple[int, ...]) -> tuple[int, int]:
    row_counts = Counter(coordinate // order for coordinate in support)
    column_counts = Counter(coordinate % order for coordinate in support)
    odd_rows = sum(count % 2 for count in row_counts.values())
    odd_columns = sum(count % 2 for count in column_counts.values())
    return odd_rows, odd_columns


def family_incidence_maxima(
    supports: tuple[tuple[int, ...], ...],
) -> tuple[int, ...]:
    values = []
    for fixed_size in range(4):
        counts: Counter[tuple[int, ...]] = Counter()
        for support in supports:
            counts.update(combinations(support, fixed_size))
        values.append(max(counts.values()))
    return tuple(values)


def exact_family_check(order: int) -> int:
    supports = tuple(combinations(range(order * order), 3))
    records = {support: support_records(order, support) for support in supports}
    families = {
        "endpoint_one": tuple(support for support in supports if records[support][1] == 1),
        "endpoint_three": tuple(support for support in supports if records[support][1] == 3),
        "middle_one_one": tuple(support for support in supports if records[support] == (1, 1)),
        "middle_one_three": tuple(
            support for support in supports if records[support] == (1, 3)
        ),
        "middle_three_one": tuple(
            support for support in supports if records[support] == (3, 1)
        ),
        "middle_three_three": tuple(support for support in supports if records[support] == (3, 3)),
    }
    expected = {
        "endpoint_one": endpoint_record_one_incidence(order),
        "endpoint_three": endpoint_record_three_incidence(order),
        "middle_one_one": middle_record_one_one_incidence(order),
        "middle_one_three": middle_record_one_three_incidence(order),
        "middle_three_one": middle_record_one_three_incidence(order),
        "middle_three_three": middle_record_three_three_incidence(order),
    }
    observed = {name: family_incidence_maxima(value) for name, value in families.items()}
    if observed != expected:
        raise AssertionError(("cubic family incidence", order, observed, expected))
    return len(supports)


def exact_q4_link_check() -> int:
    order = 4
    supports = tuple(combinations(range(order * order), 3))
    row_record = {
        support: support_records(order, support)[0] for support in supports
    }
    column_record = {
        support: support_records(order, support)[1] for support in supports
    }
    comparisons = 0
    for record, expected in (
        (1, record_one_link_bound(order)),
        (3, record_three_link_bound(order)),
    ):
        left = tuple(support for support in supports if column_record[support] == record)
        right = tuple(support for support in supports if row_record[support] == record)
        maximum = Fraction(0)
        for left_support in left:
            for right_support in right:
                value = abs(moment(order, left_support, right_support))
                maximum = max(maximum, value)
                comparisons += 1
        if maximum != expected:
            raise AssertionError(("q4 cubic link maximum", record, maximum, expected))
    return comparisons


def exact_q8_extremizer_check() -> int:
    cases = (
        ((0, 1, 8), (0, 2, 16), record_one_link_bound(8)),
        ((0, 1, 2), (0, 8, 16), record_three_link_bound(8)),
    )
    for left, right, expected in cases:
        if abs(moment(8, left, right)) != expected:
            raise AssertionError(("q8 cubic link extremizer", left, right, expected))
    return len(cases)


def main() -> None:
    candidates = candidate_entries()
    repaired = repaired_entries()
    if len(candidates) != 38 or len(repaired) != 38:
        raise AssertionError(("four-cubic repair inventory", len(candidates), len(repaired)))
    repaired_set = set(repaired)
    if any(not set(orbit(entry)).issubset(repaired_set) for entry in repaired):
        raise AssertionError("four-cubic incidence repair is not orbit closed")
    if set(coefficient_map()) != repaired_set:
        raise AssertionError("four-cubic coefficient map")
    for profile, split in repaired:
        for records in ((a, b, c) for a in (1, 3) for b in (1, 3) for c in (1, 3)):
            exact = sector_squared_coefficient(split, records)
            if Fraction.from_float(outward_sqrt(exact)) ** 2 < exact:
                raise AssertionError(("sector coefficient not outward", split, records))
        if coefficient((profile, split)) > 1:
            raise AssertionError(("four-cubic coefficient above one", split))

    q4_support_count = exact_family_check(4)
    q8_support_count = exact_family_check(8)
    link_comparisons = exact_q4_link_check()
    q8_extremizers = exact_q8_extremizer_check()
    result = diagnostic()
    if (
        result.multicubic_entries,
        result.double_cubic_entries,
        result.record_sectors,
        result.maximum_split,
        result.remaining_quarantined_entries,
    ) != (14, 24, 8, (1, 3, 0, 2), 64):
        raise AssertionError(("four-cubic diagnostic", result))
    if result.maximum_coefficient != coefficient(((3, 3, 3, 3), (1, 3, 0, 2))):
        raise AssertionError("four-cubic maximum coefficient")
    committed = (
        ROOT / "artifacts" / "q64_masked_four_cubic_incidence_repair.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale four-cubic incidence artifact")
    print(
        "q64 masked four-cubic incidence regression passed: "
        f"repaired={result.repaired_entries},"
        f"maximum={result.maximum_coefficient:.12g},"
        f"q4_supports={q4_support_count},"
        f"q8_supports={q8_support_count},"
        f"q4_link_comparisons={link_comparisons},"
        f"q8_extremizers={q8_extremizers},"
        f"remaining={result.remaining_quarantined_entries}"
    )


if __name__ == "__main__":
    main()
