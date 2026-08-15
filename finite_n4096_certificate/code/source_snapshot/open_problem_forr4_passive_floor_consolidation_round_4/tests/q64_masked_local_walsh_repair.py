#!/usr/bin/env python3
"""Regression for the 180-entry masked local-Walsh repair."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_degree_ten_completion_row_insertion import orbit  # noqa: E402
from q64_masked_local_walsh_repair import (  # noqa: E402
    artifact_text,
    candidate_entries,
    ceil_sqrt,
    coefficient,
    coefficient_map,
    diagnostic,
    entry_integer_mask_factor,
    integer_cross_mask_factor,
    mechanism,
    repaired_entries,
)
from signed_permutation_link_moment import link_moment  # noqa: E402


def walsh(order: int, left: int, right: int) -> Fraction:
    sign = -1 if (left & right).bit_count() % 2 else 1
    return Fraction(sign, order)


def endpoint_character_checks(order: int = 4) -> int:
    checks = 0
    singletons = tuple((coordinate,) for coordinate in range(order**2))
    for degree in (3, 5, 7, 9):
        for support in combinations(range(order**2), degree):
            support_xor = 0
            for coordinate in support:
                support_xor ^= coordinate
            left_weight = link_moment(order, support, singletons[0]) / walsh(
                order, support_xor, 0
            )
            right_weight = link_moment(order, singletons[0], support) / walsh(
                order, support_xor, 0
            )
            for singleton in singletons:
                expected_character = walsh(order, support_xor, singleton[0])
                if link_moment(order, support, singleton) != left_weight * expected_character:
                    raise AssertionError(("left endpoint character", degree, support, singleton))
                if link_moment(order, singleton, support) != right_weight * expected_character:
                    raise AssertionError(("right endpoint character", degree, support, singleton))
                checks += 2
    return checks


def main() -> None:
    candidates = candidate_entries()
    repaired = repaired_entries()
    if len(candidates) != 180 or candidates != repaired:
        raise AssertionError(("local Walsh repair inventory", len(candidates), len(repaired)))
    repaired_set = set(repaired)
    if any(not set(orbit(entry)).issubset(repaired_set) for entry in repaired):
        raise AssertionError("local Walsh repair is not complement/reversal closed")
    mechanisms = Counter(mechanism(entry) for entry in repaired)
    if mechanisms != {
        "internal_singleton": 80,
        "same_side_singleton_pair": 52,
        "singleton_pair_chain": 48,
    }:
        raise AssertionError(("local Walsh mechanisms", mechanisms))
    classes = Counter(
        "cubic" if entry[0] == (3, 3, 3, 3) else tuple(sorted(entry[0]))
        for entry in repaired
    )
    if classes["cubic"]:
        raise AssertionError("cubic-only entry entered local Walsh theorem")
    if (ceil_sqrt(2), integer_cross_mask_factor(1, 2), integer_cross_mask_factor(2, 3)) != (2, 3, 6):
        raise AssertionError("integer mask factor primitives")
    values = {coefficient(entry) for entry in repaired}
    if min(values) != Fraction(3, 64) or max(values) != Fraction(54, 64):
        raise AssertionError(("local Walsh coefficient range", min(values), max(values)))
    if any(coefficient(entry) != Fraction(entry_integer_mask_factor(entry), 64) for entry in repaired):
        raise AssertionError("local Walsh coefficient construction")
    if set(coefficient_map()) != repaired_set:
        raise AssertionError("local Walsh coefficient map inventory")

    exact_checks = endpoint_character_checks()
    if exact_checks != 889856:
        raise AssertionError(("q4 endpoint-character inventory", exact_checks))

    result = diagnostic()
    committed = (
        ROOT / "artifacts" / "q64_masked_local_walsh_repair.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale masked local-Walsh artifact")
    print(
        "q64 masked local-Walsh regression passed: "
        f"repaired={result.repaired_entries},"
        f"remaining={result.remaining_quarantined_entries},"
        f"maximum={result.maximum_coefficient:.12g},"
        f"q4_endpoint_checks={exact_checks}"
    )


if __name__ == "__main__":
    main()
