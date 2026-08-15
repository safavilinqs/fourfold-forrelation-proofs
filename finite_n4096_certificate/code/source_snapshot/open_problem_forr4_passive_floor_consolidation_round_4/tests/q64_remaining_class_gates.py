#!/usr/bin/env python3
"""Regression for the q64 remaining-class inventory and lead gate."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_remaining_class_gates import artifact_text, diagnostic  # noqa: E402


def main() -> None:
    result = diagnostic()
    observed_discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.balanced_open_entries,
        result.proved_entries,
        result.remaining_entries,
        result.lead_class,
        result.lead_entries,
        tuple(row.entries for row in result.classes),
    )
    expected_discrete = (
        64,
        4096,
        16_384,
        888,
        110,
        778,
        "one_split_cubic_one_split_higher",
        280,
        (280, 176, 140, 96, 48, 24, 8, 6),
    )
    if observed_discrete != expected_discrete:
        raise AssertionError(("remaining-class discrete result", observed_discrete))

    observed = (
        result.baseline.total,
        result.baseline.margin_to_one_third,
        result.lead_threshold_coefficient,
        result.lead_reserve_coefficient,
        result.lead_reserve_multiplier_over_target,
    )
    expected = (
        0.2960908671821436,
        0.0372424661511897,
        0.2255367435661101,
        0.2229211469514268,
        1.7972407795284706,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-10):
        raise AssertionError(("remaining-class numeric result", observed))

    if sum(row.entries for row in result.classes) != result.remaining_entries:
        raise AssertionError("remaining classes do not exhaust the open entries")
    if result.lead_reserve_coefficient >= result.lead_threshold_coefficient:
        raise AssertionError("reserve gate must be below threshold gate")

    committed = (ROOT / "artifacts" / "q64_remaining_class_gates.json").read_text(
        encoding="utf-8"
    )
    if committed != artifact_text(result):
        raise AssertionError("stale q64 remaining-class gate artifact")

    print(
        "q64 remaining class gates passed: "
        f"proved={result.proved_entries},"
        f"remaining={result.remaining_entries},"
        f"lead_entries={result.lead_entries},"
        f"lead_reserve={result.lead_reserve_coefficient:.12g},"
        f"lead_multiplier={result.lead_reserve_multiplier_over_target:.12g},"
        "status=routing_gate_not_theorem"
    )


if __name__ == "__main__":
    main()
