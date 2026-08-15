#!/usr/bin/env python3
"""Regression for the actual-mask quintic slice repair."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import isclose
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent
    / "open_problem_forr4_passive_floor_consolidation_round_3"
    / "searches"
)
sys.path.insert(0, str(ROOT / "searches"))
sys.path.insert(0, str(ROUND3_SEARCHES))

import occupation_compatible_sector_optimization as occupation  # noqa: E402
from q64_degree_ten_completion_row_insertion import orbit  # noqa: E402
from q64_masked_quintic_slice_repair import (  # noqa: E402
    artifact_text,
    candidate_entries,
    coefficient_map,
    complement,
    cross_distinctness_factor_upper,
    diagnostic,
    exact_quintic_singleton_slice_energies,
    repaired_entries,
    quintic_hadamard_chain_squared_bound,
    row_energy_bound,
    separated_tensor_squared_bound,
    squared_coefficient,
)
from signed_permutation_link_moment import link_moment  # noqa: E402


def record_size(support: tuple[int, ...], side: str, order: int) -> int:
    parity = [0] * order
    for coordinate in support:
        row, column = divmod(coordinate, order)
        parity[column if side == "left" else row] ^= 1
    return sum(parity)


def exact_q4_bi_record_tail_check() -> tuple[int, Fraction]:
    order = 4
    singletons = tuple((coordinate,) for coordinate in range(order**2))
    supports = tuple(
        support
        for support in combinations(range(order**2), 5)
        if record_size(support, "left", order) == 1
        and record_size(support, "right", order) == 1
    )
    tails = tuple(
        sum(link_moment(order, support, singleton) ** 2 for singleton in singletons)
        for support in supports
    )
    if set(tails) != {Fraction(1, 9)}:
        raise AssertionError(("q4 bi-record quintic tail", set(tails)))
    return len(supports), max(tails)


def walsh(order: int, left: int, right: int) -> Fraction:
    sign = -1 if (left & right).bit_count() % 2 else 1
    return Fraction(sign, order)


def quintic_endpoint_weight(
    order: int,
    support: tuple[int, ...],
    transpose: bool = False,
) -> Fraction:
    groups: dict[int, list[int]] = {}
    for coordinate in support:
        row, column = divmod(coordinate, order)
        key, value = (row, column) if transpose else (column, row)
        groups.setdefault(key, []).append(value)
    odd = [key for key, values in groups.items() if len(values) % 2]
    if len(odd) != 1:
        return Fraction(0)
    even_xors = []
    for key, values in groups.items():
        if key == odd[0]:
            continue
        value = 0
        for item in values:
            value ^= item
        even_xors.append(value)
    if not even_xors:
        return Fraction(1)
    if len(even_xors) == 1:
        return Fraction(1) if even_xors[0] == 0 else Fraction(-1, order - 1)
    if even_xors[0] == even_xors[1]:
        return Fraction(-1, order - 1)
    return Fraction(2, (order - 1) * (order - 2))


def exact_q4_quintic_endpoint_identity() -> int:
    order = 4
    supports = tuple(combinations(range(order**2), 5))
    singletons = tuple((coordinate,) for coordinate in range(order**2))
    for support in supports:
        support_xor = 0
        for coordinate in support:
            support_xor ^= coordinate
        left_weight = quintic_endpoint_weight(order, support)
        right_weight = quintic_endpoint_weight(order, support, transpose=True)
        for singleton in singletons:
            coordinate = singleton[0]
            expected = walsh(order, support_xor, coordinate)
            if link_moment(order, support, singleton) != left_weight * expected:
                raise AssertionError(("left endpoint identity", support, singleton))
            if link_moment(order, singleton, support) != right_weight * expected:
                raise AssertionError(("right endpoint identity", support, singleton))
    return len(supports) * len(singletons) * 2


def main() -> None:
    exact = exact_quintic_singleton_slice_energies(64)
    inherited = occupation.endpoint_quintic_singleton_slice_energies(64)
    if any(
        not isclose(float(left), right, rel_tol=0, abs_tol=2e-12)
        for left, right in zip(exact, inherited, strict=True)
    ):
        raise AssertionError(
            ("exact/inherited quintic slices differ", exact, inherited)
        )

    candidates = candidate_entries()
    repaired = repaired_entries()
    if len(candidates) != 54 or len(repaired) != 54:
        raise AssertionError(
            ("masked quintic inventory", len(candidates), len(repaired))
        )
    repaired_set = set(repaired)
    for entry in repaired:
        profile, split = entry
        if (profile, complement(profile, split)) not in repaired_set:
            raise AssertionError(("repair not complement closed", entry))
        if not set(orbit(entry)).issubset(repaired_set):
            raise AssertionError(("repair not orbit closed", entry))
        squared = squared_coefficient(entry)
        direct = min(
            row_energy_bound(profile, split),
            row_energy_bound(profile, complement(profile, split)),
        )
        if squared > direct or squared >= 1:
            raise AssertionError(
                ("invalid row-energy coefficient", entry, squared, direct)
            )

    expected_squared = {
        Fraction(1023, 4194304),
        Fraction(55638713045, 3670705963008),
        Fraction(1432003925, 5549064192),
        Fraction(4199425, 16257024),
        Fraction(262577446975, 721554505728),
        Fraction(791806327225, 1894080577536),
        Fraction(81, 4096),
        Fraction(81, 256),
    }
    actual_squared = {squared_coefficient(entry) for entry in repaired}
    if actual_squared != expected_squared:
        raise AssertionError(("unexpected distinct coefficients", actual_squared))

    coefficients = coefficient_map()
    if set(coefficients) != repaired_set:
        raise AssertionError("coefficient map does not equal repaired inventory")
    if any(
        Fraction.from_float(value) ** 2 < squared_coefficient(entry)
        for entry, value in coefficients.items()
    ):
        raise AssertionError("display coefficient was not rounded upward")

    separated_chain = {
        entry
        for entry in candidates
        if entry[0] == (5, 1, 1, 5)
        and entry[1][1] != entry[1][2]
        and quintic_hadamard_chain_squared_bound(entry[1])
        == squared_coefficient(entry)
    }
    if len(separated_chain) != 8:
        raise AssertionError(("unexpected separated chain", separated_chain))
    if (cross_distinctness_factor_upper(1), cross_distinctness_factor_upper(2)) != (3, 6):
        raise AssertionError("quintic mask factors")
    for entry in candidates:
        if entry[0] == (5, 1, 1, 5):
            if squared_coefficient(entry) > separated_tensor_squared_bound(entry[1]):
                raise AssertionError(("tensor bound was not retained", entry))

    q4_supports, q4_tail = exact_q4_bi_record_tail_check()
    if (q4_supports, q4_tail) != (1008, Fraction(1, 9)):
        raise AssertionError(("q4 tail inventory", q4_supports, q4_tail))
    q4_endpoint_checks = exact_q4_quintic_endpoint_identity()
    if q4_endpoint_checks != 139776:
        raise AssertionError(("q4 endpoint checks", q4_endpoint_checks))

    result = diagnostic()
    artifact = ROOT / "artifacts" / "q64_masked_quintic_slice_repair.json"
    if artifact.read_text(encoding="utf-8") != artifact_text(result):
        raise AssertionError("masked quintic artifact is stale")
    print(
        "q64 masked quintic slice repair passed: "
        f"candidates={result.candidate_entries},"
        f"repaired={result.repaired_entries},"
        f"remaining={result.remaining_quarantined_entries},"
        f"maximum_squared={result.maximum_squared_coefficient:.12g},"
        f"maximum_coefficient={result.maximum_coefficient:.12g},"
        f"q4_bi_record_tail={q4_tail},"
        f"q4_endpoint_checks={q4_endpoint_checks}"
    )


if __name__ == "__main__":
    main()
