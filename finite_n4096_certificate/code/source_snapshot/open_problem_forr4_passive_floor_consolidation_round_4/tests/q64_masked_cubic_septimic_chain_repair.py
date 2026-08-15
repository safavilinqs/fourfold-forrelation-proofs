#!/usr/bin/env python3
"""Regression for the complete masked cubic--septimic chain repair."""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_degree_ten_completion_row_insertion import orbit  # noqa: E402
from q64_masked_cubic_septimic_chain_repair import (  # noqa: E402
    artifact_text,
    candidate_entries,
    coefficient,
    coefficient_map,
    coefficient_with_mechanism,
    cubic_fixed_one_record_energies,
    degree_seven_bidegree_incidence_table,
    degree_seven_endpoint_energy_bound,
    diagnostic,
    nonzero_active_character_square_bound,
    repaired_entries,
    zero_active_completion_bound,
)
from signed_permutation_link_moment import link_moment  # noqa: E402


def record_size(support: tuple[int, ...], axis: int, order: int) -> int:
    parity = [0] * order
    for coordinate in support:
        row, column = divmod(coordinate, order)
        parity[(row, column)[axis]] ^= 1
    return sum(parity)


def direct_q4_degree_seven_incidence(
    selected: int,
    row_record: int,
    column_record: int,
) -> int:
    order = 4
    counts: dict[tuple[int, ...], int] = defaultdict(int)
    for support in combinations(range(order * order), 7):
        if (
            record_size(support, 0, order) != row_record
            or record_size(support, 1, order) != column_record
        ):
            continue
        for partial in combinations(support, selected):
            counts[partial] += 1
    return max(counts.values())


def direct_q4_cubic_fixed_one_record_energies() -> set[tuple[Fraction, Fraction]]:
    order = 4
    cells = tuple(range(order * order))
    observed = set()
    for fixed in cells:
        available = tuple(cell for cell in cells if cell != fixed)
        for singleton in cells:
            energies = {1: Fraction(0), 3: Fraction(0)}
            for completion in combinations(available, 2):
                support = tuple(sorted((fixed,) + completion))
                moment = link_moment(order, (singleton,), support)
                energies[record_size(support, 1, order)] += moment * moment
            observed.add((energies[1], energies[3]))
    return observed


def zero_active_supports(order: int) -> tuple[tuple[int, ...], ...]:
    """Generate the complete record-one/no-active septimic family."""

    zero_four_sets = tuple(
        rows
        for rows in combinations(range(order), 4)
        if rows[0] ^ rows[1] ^ rows[2] ^ rows[3] == 0
    )
    supports = set()
    for degree_three_column in range(order):
        for degree_four_column in range(order):
            if degree_three_column == degree_four_column:
                continue
            for degree_four_rows in zero_four_sets:
                for omitted in degree_four_rows:
                    degree_three_rows = tuple(
                        row for row in degree_four_rows if row != omitted
                    )
                    support = tuple(
                        sorted(
                            tuple(
                                row * order + degree_three_column
                                for row in degree_three_rows
                            )
                            + tuple(
                                row * order + degree_four_column
                                for row in degree_four_rows
                            )
                        )
                    )
                    supports.add(support)
    return tuple(sorted(supports))


def exact_zero_active_incidence(order: int, selected: int) -> int:
    counts: Counter[tuple[int, ...]] = Counter()
    for support in zero_active_supports(order):
        counts.update(combinations(support, selected))
    return max(counts.values())


def main() -> None:
    candidates = candidate_entries()
    repaired = repaired_entries()
    if len(candidates) != 12 or len(repaired) != 12:
        raise AssertionError(("cubic-septimic inventory", len(candidates), len(repaired)))
    repaired_set = set(repaired)
    if any(not set(orbit(entry)).issubset(repaired_set) for entry in repaired):
        raise AssertionError("cubic-septimic repair is not orbit closed")
    if len({frozenset(orbit(entry)) for entry in repaired}) != 3:
        raise AssertionError("cubic-septimic orbit count")

    expected_q64 = {
        (1, 3): 110_311_919,
        (1, 4): 16_209_900,
        (3, 3): 5_152_854_224,
        (3, 4): 505_158_144,
    }
    for (record, selected), expected in expected_q64.items():
        observed = degree_seven_bidegree_incidence_table(64, record, 1)[selected]
        if observed != expected:
            raise AssertionError(("q64 degree-seven incidence", record, selected, observed))

    for record in (1, 3):
        formula = degree_seven_bidegree_incidence_table(4, record, 1)
        for selected in (3, 4):
            direct = direct_q4_degree_seven_incidence(selected, record, 1)
            if formula[selected] != direct:
                raise AssertionError(
                    ("q4 degree-seven incidence", record, selected, formula[selected], direct)
                )

    fixed_one = cubic_fixed_one_record_energies(4)
    direct_energies = direct_q4_cubic_fixed_one_record_energies()
    if fixed_one != {1: Fraction(3, 16), 3: Fraction(3, 8)}:
        raise AssertionError(("q4 fixed-one formula", fixed_one))
    if direct_energies != {(fixed_one[1], fixed_one[3])}:
        raise AssertionError(("q4 fixed-one direct", direct_energies))

    zero_active_maxima = {}
    for order, expected in (
        (4, {3: 15, 4: 12}),
        (8, {3: 35, 4: 28}),
    ):
        for selected in (3, 4):
            observed = exact_zero_active_incidence(order, selected)
            zero_active_maxima[(order, selected)] = observed
            if observed != expected[selected]:
                raise AssertionError(
                    ("zero-active incidence", order, selected, observed)
                )
            if observed > zero_active_completion_bound(order):
                raise AssertionError(
                    ("zero-active analytic bound", order, selected, observed)
                )
    if nonzero_active_character_square_bound(8) != Fraction(1, 49):
        raise AssertionError("q8 nonzero-active character square")
    if nonzero_active_character_square_bound(64) != Fraction(1, 3969):
        raise AssertionError("q64 nonzero-active character square")
    if not all(
        degree_seven_endpoint_energy_bound(selected, 1)
        < degree_seven_bidegree_incidence_table(64, 1, 1)[selected]
        for selected in (3, 4)
    ):
        raise AssertionError("record-one endpoint refinement")

    coefficients = coefficient_map()
    if set(coefficients) != repaired_set or max(coefficients.values()) >= 1:
        raise AssertionError("cubic-septimic coefficient map")
    expected_coefficients = (
        0.08842196379947868,
        0.40931062674303853,
        0.5386265463506696,
    )
    observed_coefficients = tuple(sorted(set(coefficients.values())))
    if any(
        abs(left - right) > 4e-15
        for left, right in zip(observed_coefficients, expected_coefficients, strict=True)
    ):
        raise AssertionError(("cubic-septimic coefficients", observed_coefficients))
    mechanisms = Counter(
        coefficient_with_mechanism(entry)[1] for entry in repaired
    )
    if mechanisms != {
        "column_complete_cubic_endpoint": 2,
        "column_fixed_one_record_split": 4,
        "row_complete_cubic_endpoint": 2,
        "row_fixed_one_record_split": 4,
    }:
        raise AssertionError(("cubic-septimic mechanisms", mechanisms))

    result = diagnostic()
    if result.remaining_quarantined_entries != 52:
        raise AssertionError(("remaining quarantine", result.remaining_quarantined_entries))
    committed = (
        ROOT / "artifacts" / "q64_masked_cubic_septimic_chain_repair.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale cubic-septimic chain artifact")
    print(
        "q64 masked cubic-septimic chain regression passed: "
        f"repaired={result.repaired_entries},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"remaining={result.remaining_quarantined_entries},"
        "q4_incidence_checks=4,q4_endpoint_rows=256,"
        "zero_active_q4_q8=full_exact"
    )


if __name__ == "__main__":
    main()
