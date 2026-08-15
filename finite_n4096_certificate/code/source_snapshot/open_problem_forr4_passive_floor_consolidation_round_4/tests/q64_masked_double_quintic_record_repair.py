#!/usr/bin/env python3
"""Regression for the record-resolved masked double-quintic repair."""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_masked_double_quintic_record_repair import (  # noqa: E402
    Q64_FIXED_SINGLETON_FIXED_ONE,
    Q64_FIXED_SINGLETON_FIXED_TRIPLE,
    Q64_PAIR_MAXIMIZERS,
    Q64_TRIPLE_MAXIMIZERS,
    Q64_VARIABLE_SINGLETON_FIXED_PAIR,
    Q64_VARIABLE_SINGLETON_FIXED_FOUR_BOUNDS,
    artifact_text,
    diagnostic,
    enumerated_record_maxima,
    fixed_one_record_energies,
    outward_coefficient,
    repaired_entries,
    squared_coefficient_bound,
    variable_four_record_bounds,
)
from q64_masked_double_quintic_endpoint_repair import (  # noqa: E402
    candidate_entries,
    repaired_entries as endpoint_repaired_entries,
)
from signed_permutation_link_moment import moment  # noqa: E402


def records(order: int, support: tuple[int, ...]) -> tuple[int, int]:
    rows = Counter(cell // order for cell in support)
    columns = Counter(cell % order for cell in support)
    return (
        sum(value % 2 for value in rows.values()),
        sum(value % 2 for value in columns.values()),
    )


def direct_q4_maxima(
    selected: int, variable_singleton: bool
) -> dict[int, Fraction]:
    """Build every q4 endpoint row directly from physical link moments."""

    order = 4
    cells = tuple(range(order * order))
    accumulators: dict[tuple, dict[int, Fraction]] = defaultdict(
        lambda: defaultdict(Fraction)
    )
    for support in combinations(cells, 5):
        row_record, central_record = records(order, support)
        if row_record != 1:
            continue
        for singleton in cells:
            value = moment(order, (singleton,), support)
            if not value:
                continue
            for partial in combinations(support, selected):
                key = partial if variable_singleton else (partial, singleton)
                accumulators[key][central_record] += value * value
    return {
        record: max(values.get(record, Fraction(0)) for values in accumulators.values())
        for record in (1, 3, 5)
    }


def direct_selected_q8(
    partial: tuple[int, ...], variable_singleton: bool
) -> dict[int, Fraction]:
    """Directly screen one q8 partial support with exact physical moments."""

    order = 8
    cells = tuple(range(order * order))
    remaining = tuple(cell for cell in cells if cell not in partial)
    totals: dict[int, Fraction] = defaultdict(Fraction)
    for completion in combinations(remaining, 5 - len(partial)):
        support = tuple(sorted(partial + completion))
        row_record, central_record = records(order, support)
        if row_record != 1:
            continue
        if variable_singleton:
            for singleton in cells:
                value = moment(order, (singleton,), support)
                totals[central_record] += value * value
        else:
            value = moment(order, (0,), support)
            totals[central_record] += value * value
    return dict(totals)


def main() -> None:
    # The finite affine-shape formulas agree with a fully independent q4
    # enumeration of physical supports and exact signed-permutation moments.
    q4_pair = {
        record: value
        for record, value, _ in enumerated_record_maxima(4, 2, True)
    }
    q4_triple = {
        record: value
        for record, value, _ in enumerated_record_maxima(4, 3, False)
    }
    if q4_pair != direct_q4_maxima(2, True):
        raise AssertionError(("q4 pair direct screen", q4_pair))
    if q4_triple != direct_q4_maxima(3, False):
        raise AssertionError(("q4 triple direct screen", q4_triple))
    q4_fixed_one = fixed_one_record_energies(4)
    if q4_fixed_one != direct_q4_maxima(1, False):
        raise AssertionError(("q4 fixed-one direct screen", q4_fixed_one))
    q4_variable_four = direct_q4_maxima(4, True)
    q4_variable_four_bounds = variable_four_record_bounds(4)
    if any(
        q4_variable_four[record] > q4_variable_four_bounds[record]
        for record in (1, 3, 5)
    ):
        raise AssertionError(("q4 variable-four bound", q4_variable_four))

    # Selected q8 representatives test all pair affine types and the two
    # triple shapes that attain a record maximum.  These sums use the actual
    # link-moment evaluator and every cross-cut distinctness constraint.
    q8_pairs = {
        (0, 1): {1: Fraction(88, 7), 3: Fraction(908, 7), 5: Fraction(960, 7)},
        (0, 8): {1: Fraction(88, 7), 3: Fraction(220, 7)},
        (0, 9): {
            1: Fraction(204, 49),
            3: Fraction(992, 49),
            5: Fraction(960, 49),
        },
    }
    for partial, expected in q8_pairs.items():
        observed = direct_selected_q8(partial, True)
        if observed != expected:
            raise AssertionError(("q8 pair physical screen", partial, observed))
    q8_triples = {
        (0, 1, 8): {1: Fraction(3, 112), 3: Fraction(129, 1568)},
        (0, 1, 2): {
            1: Fraction(3, 448),
            3: Fraction(227, 448),
            5: Fraction(73, 112),
        },
    }
    for partial, expected in q8_triples.items():
        observed = direct_selected_q8(partial, False)
        if observed != expected:
            raise AssertionError(("q8 triple physical screen", partial, observed))

    pair = enumerated_record_maxima(64, 2, True)
    triple = enumerated_record_maxima(64, 3, False)
    if {r: v for r, v, _ in pair} != Q64_VARIABLE_SINGLETON_FIXED_PAIR:
        raise AssertionError(("q64 pair energies", pair))
    if {r: shape for r, _, shape in pair} != Q64_PAIR_MAXIMIZERS:
        raise AssertionError(("q64 pair maximizers", pair))
    if {r: v for r, v, _ in triple} != Q64_FIXED_SINGLETON_FIXED_TRIPLE:
        raise AssertionError(("q64 triple energies", triple))
    if {r: shape for r, _, shape in triple} != Q64_TRIPLE_MAXIMIZERS:
        raise AssertionError(("q64 triple maximizers", triple))
    if fixed_one_record_energies(64) != Q64_FIXED_SINGLETON_FIXED_ONE:
        raise AssertionError("q64 fixed-one record energies")
    if variable_four_record_bounds(64) != Q64_VARIABLE_SINGLETON_FIXED_FOUR_BOUNDS:
        raise AssertionError("q64 variable-four record bounds")

    expected_squares = {
        Fraction(4425307239757, 65019772670067081216),
        Fraction(119455598359945, 1462969457930881990656),
        Fraction(5901977909483, 1291082725544951808),
        Fraction(73124314015495, 1860261863620608),
    }
    repaired = repaired_entries()
    observed_squares = {squared_coefficient_bound(entry) for entry in repaired}
    if observed_squares != expected_squares:
        raise AssertionError(("double-quintic exact row energies", observed_squares))
    result = diagnostic()
    if len(repaired) != 12 or result.repaired_orbits != 5:
        raise AssertionError("double-quintic record inventory")
    if result.remaining_double_quintic_entries != 0:
        raise AssertionError("double-quintic remaining inventory")
    hard = set(candidate_entries()) - set(endpoint_repaired_entries()) - set(repaired)
    if hard:
        raise AssertionError(("double-quintic hard orbit", hard))
    for entry in repaired:
        value = outward_coefficient(entry)
        if Fraction.from_float(value) ** 2 < squared_coefficient_bound(entry):
            raise AssertionError(("coefficient not outward rounded", entry))
    committed = (
        ROOT / "artifacts" / "q64_masked_double_quintic_record_repair.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale double-quintic record artifact")
    print(
        "q64 masked double-quintic record regression passed: "
        f"repaired={result.repaired_entries},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        "q4=full_exact,q8=selected_exact,q64=finite_exact_shapes"
    )


if __name__ == "__main__":
    main()
