#!/usr/bin/env python3
"""Independent physical regression for the joint shared-quintic theorem."""

from __future__ import annotations

from collections import Counter
from decimal import Decimal, getcontext
from fractions import Fraction
from json import loads
from math import comb
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
sys.path.insert(0, str(ROOT / "searches"))

from q64_joint_recovered_cubic_quintic_contraction import (  # noqa: E402
    artifact_text,
    coefficient_map,
    diagnostic,
    feasible_quintic_pattern_pairs,
    repaired_entries,
    shape_rows,
)
from q64_recovered_cubic_quintic_independent_audit import (  # noqa: E402
    DirectQ4Plant,
    direct_permutation_moment,
    support_records,
)


EXPECTED_PATTERN_PAIRS = {
    ((2, 2, 1), (2, 2, 1)),
    ((2, 2, 1), (3, 2)),
    ((3, 2), (2, 2, 1)),
    ((2, 2, 1), (2, 1, 1, 1)),
    ((2, 2, 1), (3, 1, 1)),
    ((3, 2), (2, 1, 1, 1)),
    ((4, 1), (2, 1, 1, 1)),
    ((2, 1, 1, 1), (2, 2, 1)),
    ((2, 1, 1, 1), (3, 2)),
    ((2, 1, 1, 1), (4, 1)),
    ((3, 1, 1), (2, 2, 1)),
    ((2, 1, 1, 1), (2, 1, 1, 1)),
    ((2, 1, 1, 1), (3, 1, 1)),
    ((3, 1, 1), (2, 1, 1, 1)),
    ((3, 1, 1), (3, 1, 1)),
}


def axis_pattern(support: tuple[int, ...], order: int, axis: int) -> tuple[int, ...]:
    counts = Counter(divmod(cell, order)[axis] for cell in support)
    return tuple(sorted(counts.values(), reverse=True))


def manual_record_one_bound(order: int) -> Fraction:
    return Fraction(order + 2, order * (order - 1) * (order - 2))


def manual_plain_record_three(order: int) -> Fraction:
    return Fraction(1, comb(order, 3))


def manual_even_record_three(order: int) -> Fraction:
    return min(
        manual_plain_record_three(order),
        Fraction(3, (order - 3) * comb(order, 3)),
    )


def manual_shape_bound(
    order: int,
    row_pattern: tuple[int, ...],
    column_pattern: tuple[int, ...],
) -> Fraction:
    row_record = sum(value % 2 for value in row_pattern)
    column_record = sum(value % 2 for value in column_pattern)
    if row_record == 1:
        left = (
            Fraction(1, order * (order - 1))
            * manual_record_one_bound(order)
        )
    elif row_pattern == (3, 1, 1):
        left = Fraction(1, order) * manual_plain_record_three(order)
    elif row_pattern == (2, 1, 1, 1):
        left = Fraction(1, order) * manual_even_record_three(order)
    else:
        raise AssertionError(("unsupported manual row pattern", row_pattern))
    if column_record == 1:
        right = (
            Fraction(1, order)
            if column_pattern in ((5,), (4, 1))
            else manual_record_one_bound(order)
        )
    elif column_pattern == (3, 1, 1):
        right = manual_plain_record_three(order)
    elif column_pattern == (2, 1, 1, 1):
        right = manual_even_record_three(order)
    else:
        raise AssertionError(("unsupported manual column pattern", column_pattern))
    return left * right


def q4_full_physical_joint_maxima() -> tuple[int, int]:
    """Enumerate every q4 quintic and both adjacent cubic families exactly."""

    plant = DirectQ4Plant()
    group = plant.group_size
    endpoint = np.abs(plant.moments(1, 3).astype(np.int64)).max(axis=0)
    left = np.abs(plant.moments(3, 5).astype(np.int64))
    right = np.abs(plant.moments(5, 3).astype(np.int64))
    cubics = {
        record: np.asarray(
            [
                index
                for index, support in enumerate(plant.supports[3])
                if support_records(support, 4) == (1, record)
            ],
            dtype=np.int64,
        )
        for record in (1, 3)
    }
    terminal_cubics = {
        record: np.asarray(
            [
                index
                for index, support in enumerate(plant.supports[3])
                if support_records(support, 4)[0] == record
            ],
            dtype=np.int64,
        )
        for record in (1, 3)
    }
    quintics = {
        records: np.asarray(
            [
                index
                for index, support in enumerate(plant.supports[5])
                if support_records(support, 4) == records
            ],
            dtype=np.int64,
        )
        for records in ((1, 1), (1, 3), (3, 1), (3, 3))
    }
    expected = {
        (1, 1): Fraction(1, 192),
        (1, 3): Fraction(1, 192),
        (3, 1): Fraction(1, 64),
        (3, 3): Fraction(1, 64),
    }
    observed_patterns = set()
    for support in plant.supports[5]:
        records = support_records(support, 4)
        if records in quintics:
            observed_patterns.add(
                (axis_pattern(support, 4, 0), axis_pattern(support, 4, 1))
            )
    if observed_patterns != EXPECTED_PATTERN_PAIRS:
        raise AssertionError(("q4 quintic shape inventory", observed_patterns))
    for records, q_indices in quintics.items():
        first, second = records
        c_indices = cubics[first]
        d_indices = terminal_cubics[second]
        left_products = (
            endpoint[c_indices, None]
            * left[np.ix_(c_indices, q_indices)]
        )
        left_maxima = left_products.max(axis=0)
        right_maxima = right[np.ix_(q_indices, d_indices)].max(axis=1)
        direct = Fraction(
            int((left_maxima * right_maxima).max()),
            group**3,
        )
        if direct != expected[records]:
            raise AssertionError(("q4 joint physical maximum", records, direct))
    return len(observed_patterns), len(expected)


Q8_CASES = (
    # records, row pattern, column pattern, quintic, cubic, terminal cubic,
    # endpoint moment, left middle moment, right middle moment
    ((1, 1), (2, 2, 1), (2, 2, 1), (0, 1, 8, 9, 18), (0, 1, 32), (0, 1, 2), -Fraction(1, 56), Fraction(3, 280), Fraction(1, 56)),
    ((1, 1), (2, 2, 1), (3, 2), (0, 1, 8, 9, 16), (0, 1, 32), (0, 2, 16), -Fraction(1, 56), Fraction(3, 280), -Fraction(5, 168)),
    ((1, 1), (3, 2), (2, 2, 1), (0, 1, 2, 8, 9), (0, 2, 16), (0, 1, 2), -Fraction(1, 56), -Fraction(5, 168), Fraction(1, 56)),
    ((1, 3), (2, 2, 1), (2, 1, 1, 1), (0, 1, 8, 10, 20), (0, 4, 32), (0, 19, 32), -Fraction(1, 56), Fraction(19, 840), Fraction(3, 280)),
    ((1, 3), (2, 2, 1), (3, 1, 1), (0, 1, 8, 10, 16), (0, 4, 32), (0, 9, 19), -Fraction(1, 56), -Fraction(19, 840), Fraction(1, 56)),
    ((1, 3), (3, 2), (2, 1, 1, 1), (0, 1, 2, 8, 12), (0, 1, 16), (0, 16, 32), -Fraction(1, 56), -Fraction(5, 168), -Fraction(3, 280)),
    ((1, 3), (4, 1), (2, 1, 1, 1), (0, 1, 2, 4, 8), (0, 3, 16), (0, 16, 32), Fraction(1, 56), Fraction(5, 168), -Fraction(3, 280)),
    ((3, 1), (2, 1, 1, 1), (2, 2, 1), (0, 1, 8, 17, 34), (0, 2, 4), (0, 4, 32), Fraction(1, 8), -Fraction(1, 280), Fraction(19, 840)),
    ((3, 1), (2, 1, 1, 1), (3, 2), (0, 1, 8, 16, 33), (0, 2, 4), (0, 2, 8), Fraction(1, 8), -Fraction(3, 280), -Fraction(5, 168)),
    ((3, 1), (2, 1, 1, 1), (4, 1), (0, 1, 8, 16, 32), (0, 2, 4), (0, 2, 24), Fraction(1, 8), -Fraction(3, 280), Fraction(5, 168)),
    ((3, 1), (3, 1, 1), (2, 2, 1), (0, 1, 2, 8, 17), (0, 1, 2), (0, 4, 32), Fraction(1, 8), -Fraction(1, 168), -Fraction(19, 840)),
    ((3, 3), (2, 1, 1, 1), (2, 1, 1, 1), (0, 1, 8, 18, 36), (0, 2, 4), (0, 18, 35), Fraction(1, 8), -Fraction(1, 280), Fraction(3, 280)),
    ((3, 3), (2, 1, 1, 1), (3, 1, 1), (0, 1, 8, 16, 34), (0, 1, 4), (0, 11, 19), Fraction(1, 8), -Fraction(1, 280), Fraction(1, 56)),
    ((3, 3), (3, 1, 1), (2, 1, 1, 1), (0, 1, 2, 8, 20), (0, 1, 2), (0, 19, 32), Fraction(1, 8), -Fraction(1, 168), Fraction(3, 280)),
    ((3, 3), (3, 1, 1), (3, 1, 1), (0, 1, 2, 8, 16), (0, 3, 4), (0, 11, 19), -Fraction(1, 8), Fraction(1, 56), Fraction(1, 56)),
)


def q8_direct_shape_checks() -> tuple[int, Fraction, Fraction]:
    cache: dict[tuple[tuple[int, ...], tuple[int, ...]], Fraction] = {}

    def direct(left: tuple[int, ...], right: tuple[int, ...]) -> Fraction:
        key = left, right
        if key not in cache:
            cache[key] = direct_permutation_moment(8, left, right)
        return cache[key]

    singleton = (63,)
    plus_actual = None
    for (
        records,
        row_pattern,
        column_pattern,
        quintic,
        cubic,
        terminal,
        endpoint_expected,
        left_expected,
        right_expected,
    ) in Q8_CASES:
        if support_records(quintic, 8) != records:
            raise AssertionError(("q8 record sector", quintic))
        if axis_pattern(quintic, 8, 0) != row_pattern or axis_pattern(
            quintic, 8, 1
        ) != column_pattern:
            raise AssertionError(("q8 physical shape", quintic))
        observed = (
            direct(singleton, cubic),
            direct(cubic, quintic),
            direct(quintic, terminal),
        )
        expected = endpoint_expected, left_expected, right_expected
        if observed != expected:
            raise AssertionError(("q8 direct permutation moments", quintic, observed))
        actual = abs(observed[0] * observed[1] * observed[2])
        if actual > manual_shape_bound(8, row_pattern, column_pattern):
            raise AssertionError(("q8 joint shape bound", quintic, actual))
        if row_pattern == column_pattern == (3, 1, 1):
            plus_actual = actual

    # The endpoint-compatible 41 shape must also cover a zero xor in its
    # even four-group; this is distinct from the old incompatible 1/q link.
    zero_xor_compatible = direct(
        (0, 1, 8),
        (0, 1, 2, 3, 8),
    )
    if abs(zero_xor_compatible) != Fraction(1, 56):
        raise AssertionError(
            ("q8 endpoint-compatible 41 zero xor", zero_xor_compatible)
        )
    if abs(zero_xor_compatible) > manual_record_one_bound(8):
        raise AssertionError("q8 endpoint-compatible 41 bound")

    # Retain the exact geometry that rejected the independent-maxima proof.
    counterexample = (
        direct((0,), (0, 1, 2))
        * direct((0, 1, 2), (0, 1, 8, 16, 24))
        * direct((0, 1, 8, 16, 24), (0, 1, 2))
    )
    if abs(counterexample) != Fraction(1, 17920):
        raise AssertionError(("lost rejected counterexample", counterexample))
    expected_plus = Fraction(1, 25088)
    if plus_actual != expected_plus:
        raise AssertionError(("q8 sharp plus shape", plus_actual))
    return len(Q8_CASES), abs(counterexample), expected_plus


def q64_artifact_audit() -> tuple[int, float]:
    result = diagnostic()
    generated = artifact_text(result)
    artifact = loads(generated)
    committed = (
        ROOT
        / "artifacts"
        / "q64_joint_recovered_cubic_quintic_contraction.json"
    ).read_text(encoding="utf-8")
    if committed != generated:
        raise AssertionError("stale joint contraction artifact")
    if artifact["schema"] != "round4_q64_joint_recovered_cubic_quintic_contraction_v1":
        raise AssertionError("joint contraction schema")
    if len(artifact["shape_table"]) != 15 or len(artifact["registry_entries"]) != 12:
        raise AssertionError("joint contraction artifact inventory")
    getcontext().prec = 100
    maximum = 0.0
    for row in artifact["registry_entries"]:
        exact_sum = Decimal(0)
        for sector in row["sectors"]:
            exact = (
                Decimal(sector["squared_coefficient_numerator"])
                / Decimal(sector["squared_coefficient_denominator"])
            )
            exact_sum += exact.sqrt()
        displayed = Decimal.from_float(row["outward_coefficient"])
        if displayed < exact_sum:
            raise AssertionError(("joint coefficient rounded inward", row))
        maximum = max(maximum, row["outward_coefficient"])
    expected = (
        0.11863696369021283,
        0.3144332243430522,
        0.33828697324447987,
    )
    if tuple(sorted(set(coefficient_map().values()))) != expected:
        raise AssertionError(("q64 joint coefficients", coefficient_map()))
    if len(repaired_entries()) != 12 or maximum >= 1 or 1 - maximum < 0.6:
        raise AssertionError(("q64 joint theorem margin", maximum))
    if {
        (row.row_pattern, row.column_pattern) for row in shape_rows()
    } != set(feasible_quintic_pattern_pairs()):
        raise AssertionError("production shape table mismatch")
    return len(artifact["registry_entries"]), maximum


def main() -> None:
    q4_shapes, q4_sectors = q4_full_physical_joint_maxima()
    q8_shapes, counterexample, sharp_plus = q8_direct_shape_checks()
    q64_entries, maximum = q64_artifact_audit()
    print(
        "q64 joint recovered cubic-quintic regression passed: "
        f"q4_shapes={q4_shapes},q4_sectors={q4_sectors},"
        f"q8_shape_checks={q8_shapes},"
        f"counterexample={counterexample},"
        f"sharp_plus={sharp_plus},"
        f"q64_entries={q64_entries},maximum={maximum:.12g},"
        "status=joint_shared_quintic_theorem"
    )


if __name__ == "__main__":
    main()
