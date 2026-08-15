#!/usr/bin/env python3
"""Regression for the universal q64 multicubic theorem insertion."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_remaining_class_gates import contraction_class  # noqa: E402
from q64_universal_multicubic_insertion import (  # noqa: E402
    artifact_text,
    diagnostic,
    multicubic_entries,
)


def main() -> None:
    result = diagnostic()
    observed_discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.three_split_cubic_entries,
        result.four_split_cubic_entries,
        result.newly_closed_entries,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
    )
    expected_discrete = (64, 4096, 16_384, 8, 6, 14, 206, 220, 668)
    if observed_discrete != expected_discrete:
        raise AssertionError(("universal multicubic discrete result", observed_discrete))

    observed = (
        result.universal_coefficient,
        result.previous_routing.total,
        result.multicubic_inserted.total,
        result.multicubic_inserted.beta,
        result.multicubic_inserted.perron_upper,
        result.multicubic_inserted.promise_loss,
        result.multicubic_inserted.margin_to_one_third,
        result.margin_spent,
        result.reserve_after_declared_allowance,
    )
    expected = (
        1.0,
        0.3293832216221608,
        0.33193582943438404,
        0.7460717495549545,
        0.31456527793840144,
        0.017370551495982593,
        0.0013975038989492705,
        0.0025526078122232176,
        0.0003975038989492705,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("universal multicubic numeric result", observed))

    entries = multicubic_entries()
    if len(entries) != 14 or {
        contraction_class(entry) for entry in entries
    } != {"three_split_cubics", "four_split_cubics"}:
        raise AssertionError("universal multicubic entry class")

    committed = (
        ROOT / "artifacts" / "q64_universal_multicubic_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 universal multicubic artifact")

    print(
        "q64 universal multicubic insertion passed: "
        f"entries={result.newly_closed_entries},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.multicubic_inserted.total:.12g},"
        f"margin={result.multicubic_inserted.margin_to_one_third:.12g},"
        "status=quarantined_unmasked_coefficient_one_diagnostic"
    )


if __name__ == "__main__":
    main()
