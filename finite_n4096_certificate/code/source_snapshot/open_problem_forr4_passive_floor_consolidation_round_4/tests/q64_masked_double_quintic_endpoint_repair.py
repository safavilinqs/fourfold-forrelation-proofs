#!/usr/bin/env python3
"""Regression for the six-entry double-quintic endpoint repair."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_degree_ten_completion_row_insertion import orbit  # noqa: E402
from q64_masked_double_quintic_endpoint_repair import (  # noqa: E402
    artifact_text,
    candidate_entries,
    coefficient_map,
    diagnostic,
    outward_coefficient,
    repaired_entries,
    squared_coefficient,
)


def main() -> None:
    candidates = candidate_entries()
    repaired = repaired_entries()
    if len(candidates) != 18 or len(repaired) != 6:
        raise AssertionError(("double quintic inventory", len(candidates), len(repaired)))
    repaired_set = set(repaired)
    if any(not set(orbit(entry)).issubset(repaired_set) for entry in repaired):
        raise AssertionError("double quintic endpoint repair is not orbit closed")
    values = {squared_coefficient(entry) for entry in repaired}
    expected = {
        Fraction(1023, 1024),
        Fraction(1046529, 1048576),
    }
    if values != expected or max(values) != Fraction(1023, 1024):
        raise AssertionError(("double quintic coefficients", values))
    if any(Fraction.from_float(outward_coefficient(entry)) ** 2 < squared_coefficient(entry) for entry in repaired):
        raise AssertionError("double quintic coefficient not rounded outward")
    if set(coefficient_map()) != repaired_set:
        raise AssertionError("double quintic coefficient map")

    result = diagnostic()
    if (
        result.one_whole_endpoint_entries,
        result.oriented_double_one_four_entries,
    ) != (4, 2):
        raise AssertionError("double quintic structural inventory")
    committed = (
        ROOT / "artifacts" / "q64_masked_double_quintic_endpoint_repair.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale double quintic endpoint artifact")
    print(
        "q64 masked double-quintic endpoint regression passed: "
        f"repaired={result.repaired_entries},"
        f"maximum={result.maximum_coefficient:.12g},"
        f"remaining={result.remaining_quarantined_entries}"
    )


if __name__ == "__main__":
    main()
