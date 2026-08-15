#!/usr/bin/env python3
"""Regression for the q64 universal double-cubic insertion."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_universal_double_cubic_insertion import (  # noqa: E402
    TARGET_CLASS,
    artifact_text,
    diagnostic,
    double_cubic_entries,
)
from q64_remaining_class_gates import contraction_class  # noqa: E402


def main() -> None:
    result = diagnostic()
    observed_discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.newly_closed_entries,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
    )
    expected_discrete = (64, 4096, 16_384, 24, 236, 260, 628)
    if observed_discrete != expected_discrete:
        raise AssertionError(("double-cubic discrete result", observed_discrete))

    observed = (
        result.universal_coefficient,
        result.previous_routing.total,
        result.double_cubic_inserted.total,
        result.double_cubic_inserted.beta,
        result.double_cubic_inserted.perron_upper,
        result.double_cubic_inserted.promise_loss,
        result.double_cubic_inserted.margin_to_one_third,
        result.margin_spent,
        result.reserve_after_declared_allowance,
    )
    expected = (
        1.0,
        0.32303469500415644,
        0.33056335386676317,
        0.7460779152389242,
        0.3132184264212764,
        0.01734492744548678,
        0.0027699794665701494,
        0.00752865886260673,
        0.0017699794665701494,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("double-cubic numeric result", observed))

    entries = double_cubic_entries()
    if len(entries) != 24 or any(
        contraction_class(entry) != TARGET_CLASS for entry in entries
    ):
        raise AssertionError("double-cubic entry partition")

    committed = (
        ROOT / "artifacts" / "q64_universal_double_cubic_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 universal double-cubic artifact")

    print(
        "q64 universal double-cubic insertion passed: "
        f"entries={result.newly_closed_entries},"
        f"total={result.double_cubic_inserted.total:.12g},"
        f"margin={result.double_cubic_inserted.margin_to_one_third:.12g},"
        f"reserve_after_allowance={result.reserve_after_declared_allowance:.12g},"
        "status=quarantined_unmasked_coefficient_one_diagnostic"
    )


if __name__ == "__main__":
    main()
