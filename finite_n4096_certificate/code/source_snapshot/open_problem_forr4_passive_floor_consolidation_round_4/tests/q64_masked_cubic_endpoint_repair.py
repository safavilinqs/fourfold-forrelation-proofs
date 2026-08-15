#!/usr/bin/env python3
"""Regression for the 12-entry masked cubic-endpoint repair."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_degree_ten_completion_row_insertion import orbit  # noqa: E402
from q64_masked_cubic_endpoint_repair import (  # noqa: E402
    artifact_text,
    candidate_entries,
    coefficient,
    coefficient_map,
    cubic_endpoint_squared_factor,
    diagnostic,
    entry_kind,
    repaired_entries,
    squared_coefficient,
)
from signed_permutation_link_moment import link_moment  # noqa: E402


def exact_q4_fixed_pair_energy() -> tuple[int, set[Fraction]]:
    order = 4
    cells = tuple(range(order**2))
    values = set()
    rows = 0
    for pair in combinations(cells, 2):
        for singleton in cells:
            values.add(
                sum(
                    link_moment(
                        order,
                        tuple(sorted(pair + (extra,))),
                        (singleton,),
                    )
                    ** 2
                    for extra in cells
                    if extra not in pair
                )
            )
            rows += 1
    return rows, values


def main() -> None:
    candidates = candidate_entries()
    repaired = repaired_entries()
    if len(candidates) != 12 or candidates != repaired:
        raise AssertionError(("cubic endpoint inventory", len(candidates), len(repaired)))
    repaired_set = set(repaired)
    if any(not set(orbit(entry)).issubset(repaired_set) for entry in repaired):
        raise AssertionError("cubic endpoint repair is not orbit closed")
    kinds = {kind: sum(entry_kind(entry) == kind for entry in repaired) for kind in {
        "cubic_septimic", "recovered_cubic_quintic"
    }}
    if kinds != {"cubic_septimic": 4, "recovered_cubic_quintic": 8}:
        raise AssertionError(("cubic endpoint kinds", kinds))
    exact = Fraction(225, 4) * Fraction(3970, 258048)
    if cubic_endpoint_squared_factor() != Fraction(3970, 258048):
        raise AssertionError("q64 cubic endpoint exact factor")
    if squared_coefficient() != exact or exact >= 1:
        raise AssertionError(("cubic endpoint coefficient", squared_coefficient()))
    if Fraction.from_float(coefficient()) ** 2 < exact:
        raise AssertionError("cubic endpoint coefficient not rounded outward")
    if set(coefficient_map()) != repaired_set:
        raise AssertionError("cubic endpoint coefficient map")

    rows, q4_values = exact_q4_fixed_pair_energy()
    if rows != 1920 or q4_values != {Fraction(1, 24), Fraction(5, 24)}:
        raise AssertionError(("q4 fixed-pair cubic energy", rows, q4_values))
    if max(q4_values) != cubic_endpoint_squared_factor(4):
        raise AssertionError("q4 formula/direct cubic endpoint mismatch")

    result = diagnostic()
    committed = (
        ROOT / "artifacts" / "q64_masked_cubic_endpoint_repair.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale masked cubic-endpoint artifact")
    print(
        "q64 masked cubic-endpoint regression passed: "
        f"repaired={result.repaired_entries},"
        f"coefficient={result.coefficient:.12g},"
        f"remaining={result.remaining_quarantined_entries},"
        f"q4_rows={rows}"
    )


if __name__ == "__main__":
    main()
