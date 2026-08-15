#!/usr/bin/env python3
"""Regression for the sixth balanced whole-row contraction."""

from __future__ import annotations

from itertools import combinations
from math import comb
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from adjacent_balanced_row_slice_contraction import (  # noqa: E402
    adjacent_balanced_row_coefficient,
    adjacent_balanced_row_contraction,
    exact_fixed_row_energy,
    exact_q4_row_slice,
    record_one_row_energy_bound,
    record_three_cubic_weight,
    record_three_tail_bounds,
    target_orbit_entries,
)
from adjacent_cubic_quintic_mixed_orbit_q4 import (  # noqa: E402
    combined_link_moment,
)
from adjacent_cubic_quintic_orbit_witness import (  # noqa: E402
    parity_record_size,
    record_one_link_moment,
    unnormalized_sylvester,
)


SEED = 2026071506


def row_multiplicities(order: int, support: tuple[int, ...]) -> tuple[int, ...]:
    counts: dict[int, int] = {}
    for cell in support:
        row = cell // order
        counts[row] = counts.get(row, 0) + 1
    return tuple(sorted(counts.values(), reverse=True))


def extension_count_checks(order: int) -> None:
    """Enumerate the fixed-triple row-pattern counts used by the proof."""

    q = order
    dimension = q * q
    representatives = {
        "3": (0, 1, 2),
        "21": (0, 1, q),
        "111": (0, q, 2 * q),
    }
    expected_record_one = {
        "3": comb(q - 3, 2) + (q - 1) * comb(q, 2) + (q - 3) * (q - 1) * q,
        "21": (q - 2) * (3 * q * q + q - 6) // 2,
        "111": 3 * (q - 1) ** 2,
    }
    expected_record_three = {
        "3": (comb(q - 1, 2) * q * q, 0),
        "21": (q * (q - 2) ** 2, q * q * comb(q - 2, 2)),
        "111": (
            3 * comb(q - 1, 2),
            (q - 3) * comb(q, 2) + 3 * (q - 1) * (q - 3) * q,
        ),
    }
    for key, fixed in representatives.items():
        fixed_set = set(fixed)
        available = tuple(cell for cell in range(dimension) if cell not in fixed_set)
        record_one = 0
        no_even = 0
        one_even = 0
        for pair in combinations(available, 2):
            support = tuple(sorted(fixed + pair))
            record = parity_record_size(q, support, axis=0)
            if record == 1:
                record_one += 1
            if record != 3:
                continue
            multiplicities = row_multiplicities(q, support)
            even_groups = sum(value % 2 == 0 for value in multiplicities)
            if even_groups == 0:
                no_even += 1
            elif even_groups == 1:
                one_even += 1
            else:
                raise AssertionError(("record-three row pattern", support))
        if record_one != expected_record_one[key]:
            raise AssertionError(("record-one extensions", q, key, record_one))
        if (no_even, one_even) != expected_record_three[key]:
            raise AssertionError(("record-three extensions", q, key, no_even, one_even))


def horizontal_tail_checks(order: int) -> None:
    """Check that every record-three cubic has the exact horizontal tail."""

    q = order
    dimension = q * q
    fixed_cell = 0
    fixed_triple = (0, 1, 2)
    fixed_set = set(fixed_triple)
    quintics = tuple(
        tuple(sorted(fixed_triple + pair))
        for pair in combinations(
            tuple(cell for cell in range(dimension) if cell not in fixed_set),
            2,
        )
        if parity_record_size(
            q,
            tuple(sorted(fixed_triple + pair)),
            axis=0,
        )
        == 3
    )
    cubics = tuple(
        tuple(sorted((fixed_cell,) + pair))
        for pair in combinations(range(1, dimension), 2)
        if parity_record_size(
            q,
            tuple(sorted((fixed_cell,) + pair)),
            axis=0,
        )
        == 1
        and parity_record_size(
            q,
            tuple(sorted((fixed_cell,) + pair)),
            axis=1,
        )
        == 3
    )
    expected = 3 / ((q - 1) * (q - 2))
    for cubic in cubics:
        observed = sum(
            combined_link_moment(q, cubic, quintic) ** 2 for quintic in quintics
        )
        if not np.isclose(observed, expected, atol=3e-13):
            raise AssertionError(("horizontal record-three tail", q, cubic))


def direct_sparse_tensor_checks(q4_coefficient: float) -> float:
    """Stress the complete indexed target under correlated diagonal laws."""

    rng = np.random.default_rng(SEED)
    order = 4
    dimension = order * order
    pairs = tuple(combinations(range(dimension), 2))
    triples = tuple(combinations(range(dimension), 3))
    rows = tuple((b, 0, fixed) for b in range(dimension) for fixed in triples)
    columns = tuple(
        (a, cubic_pair, quintic_pair)
        for a in range(dimension)
        for cubic_pair in pairs
        for quintic_pair in pairs
    )
    hadamard = unnormalized_sylvester(dimension).astype(float) / np.sqrt(dimension)
    worst = 0.0
    for _ in range(10):
        selected_rows = tuple(
            rows[index] for index in rng.choice(len(rows), size=22, replace=False)
        )
        selected_columns = tuple(
            columns[index] for index in rng.choice(len(columns), size=34, replace=False)
        )
        tensor = np.zeros((len(selected_rows), len(selected_columns)))
        for row_index, (b, x, fixed) in enumerate(selected_rows):
            fixed_set = set(fixed)
            for column_index, (a, cubic_pair, quintic_pair) in enumerate(
                selected_columns
            ):
                if x in cubic_pair or fixed_set.intersection(quintic_pair):
                    continue
                cubic = tuple(sorted((x,) + cubic_pair))
                quintic = tuple(sorted(fixed + quintic_pair))
                tensor[row_index, column_index] = (
                    hadamard[b, a]
                    * record_one_link_moment(order, (b,), cubic)
                    * combined_link_moment(order, cubic, quintic)
                )
        row_law = rng.dirichlet(np.ones(len(selected_rows)))
        column_law = rng.dirichlet(np.ones(len(selected_columns)))
        weighted = np.sqrt(row_law)[:, None] * tensor * np.sqrt(column_law)[None, :]
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > q4_coefficient * (1 + 5e-12):
            raise AssertionError(("sparse target tensor", nuclear, q4_coefficient))
        worst = max(worst, nuclear / q4_coefficient)
    return worst


def main() -> None:
    extension_count_checks(4)
    extension_count_checks(8)
    horizontal_tail_checks(4)
    horizontal_tail_checks(8)

    exact_q4 = exact_q4_row_slice()
    if not np.isclose(
        exact_q4.maximum_row_energy,
        0.3763020833333333,
        atol=2e-14,
    ):
        raise AssertionError(("q4 maximum row energy", exact_q4))
    if exact_q4.maximizing_triple != (4, 8, 12):
        raise AssertionError(("q4 maximizing triple", exact_q4))

    q8_horizontal = exact_fixed_row_energy(8, 0, 0, (0, 1, 2))
    if not np.isclose(
        q8_horizontal.record_one_energy,
        0.0093143563228863,
        atol=3e-14,
    ):
        raise AssertionError(("q8 horizontal record one", q8_horizontal))
    if not np.isclose(
        q8_horizontal.record_three_energy,
        0.0334821428571429,
        atol=3e-14,
    ):
        raise AssertionError(("q8 horizontal record three", q8_horizontal))
    if not record_one_row_energy_bound(8) > q8_horizontal.record_one_energy:
        raise AssertionError(("record-one theorem bound", q8_horizontal))
    if not np.isclose(
        record_three_cubic_weight(8) * record_three_tail_bounds(8)[0],
        q8_horizontal.record_three_energy,
        atol=3e-14,
    ):
        raise AssertionError(("record-three horizontal formula", q8_horizontal))

    orbit = target_orbit_entries()
    if len(orbit) != 4:
        raise AssertionError(("target orbit", orbit))
    coefficient = adjacent_balanced_row_coefficient()
    if not np.isclose(coefficient, 0.0422410016249075, atol=3e-15):
        raise AssertionError(("q32 coefficient", coefficient))
    result = adjacent_balanced_row_contraction()
    if not result.coefficient > result.provisional_coefficient:
        raise AssertionError(("safe coefficient should exceed provisional", result))
    if not result.coefficient < result.acceptance_gate:
        raise AssertionError(("sixth-orbit gate", result))
    if not np.isclose(
        result.optimized_total,
        0.332775779206186,
        atol=4e-11,
    ):
        raise AssertionError(("six-theorem ledger", result))
    if not result.threshold_slack > 0.000557:
        raise AssertionError(("six-theorem slack", result))
    if result.next_unresolved_entries[0] != (
        (3, 1, 1, 5),
        (0, 1, 1, 3),
    ):
        raise AssertionError(("reranked next orbit", result))
    if not np.isclose(
        result.next_admissible_coefficient,
        0.0484819899411186,
        atol=4e-14,
    ):
        raise AssertionError(("next coefficient gate", result))

    worst = direct_sparse_tensor_checks(exact_q4.maximum_coefficient)
    print(
        "adjacent balanced row-slice contraction passed: "
        f"q4_energy={exact_q4.maximum_row_energy:.12g},"
        f"q8_horizontal={q8_horizontal.total_energy:.12g},"
        f"q32_record_one={result.record_one_row_energy_bound:.12g},"
        f"q32_record_three={result.record_three_row_energy_bound:.12g},"
        f"coefficient={result.coefficient:.12g},"
        f"ledger_total={result.optimized_total:.12g},"
        f"threshold_slack={result.threshold_slack:.12g},"
        f"next_admissible={result.next_admissible_coefficient:.12g},"
        f"q4_sparse_ratio={worst:.12g}"
    )


if __name__ == "__main__":
    main()
