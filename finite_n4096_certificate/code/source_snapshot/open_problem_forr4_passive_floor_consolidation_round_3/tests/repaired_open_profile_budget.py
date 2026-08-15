#!/usr/bin/env python3
"""Regression for the repaired high-degree profile budget."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from repaired_open_profile_budget import (
    ADJACENT_SPLIT_ENTRIES,
    CERTIFIED_LEADING_COEFFICIENT,
    LEADING_SPLIT_ENTRIES,
    WITNESS_COEFFICIENT,
    adjacent_forced_coefficients,
    coarse_open_completion_coefficients,
    coarse_open_completion_target,
    extended_open_budget,
    forced_coefficients,
    hybrid_open_budget,
    open_profiles,
    repaired_open_budget,
    reversal_orbits,
)


def main() -> None:
    profiles = open_profiles()
    if len(profiles) != 43 or len(reversal_orbits(profiles)) != 23:
        raise AssertionError(("open profile inventory", len(profiles)))

    repaired = repaired_open_budget()
    if not np.isclose(repaired.baseline_total, 0.283025649927006, atol=3e-12):
        raise AssertionError(("repaired baseline", repaired))
    if not np.isclose(
        repaired.common_coefficient_threshold,
        0.0367450372915,
        atol=3e-12,
    ):
        raise AssertionError(("repaired common threshold", repaired))

    hybrid = hybrid_open_budget()
    if not np.isclose(hybrid.baseline_total, 0.283590836856527, atol=3e-12):
        raise AssertionError(("hybrid baseline", hybrid))
    if not np.isclose(
        hybrid.common_coefficient_threshold,
        0.0382742303294,
        atol=3e-12,
    ):
        raise AssertionError(("hybrid common threshold", hybrid))

    extended = extended_open_budget()
    if not np.isclose(extended.baseline_total, 0.283376978283256, atol=3e-12):
        raise AssertionError(("extended baseline", extended))
    if not np.isclose(
        extended.common_coefficient_threshold,
        0.0383748987121,
        atol=3e-12,
    ):
        raise AssertionError(("extended common threshold", extended))

    # The three physical orbits must survive the coarse mapping unchanged.
    first = forced_coefficients()
    forced = adjacent_forced_coefficients()
    if len(first) != 4 or len(forced) != 12:
        raise AssertionError(("forced witness count", len(first), len(forced)))
    for entry in first:
        if not np.isclose(forced[entry], WITNESS_COEFFICIENT, atol=1e-15):
            raise AssertionError(("first witness changed", entry, forced[entry]))
    for entry in LEADING_SPLIT_ENTRIES:
        if not np.isclose(
            forced[entry], CERTIFIED_LEADING_COEFFICIENT, atol=1e-15
        ):
            raise AssertionError(("second witness changed", entry, forced[entry]))
    if any(entry not in forced for entry in ADJACENT_SPLIT_ENTRIES):
        raise AssertionError("adjacent safe orbit missing")

    coarse_coefficients = coarse_open_completion_coefficients()
    if len(coarse_coefficients) != 6016:
        raise AssertionError(("coarse coefficient inventory", len(coarse_coefficients)))
    if any(coarse_coefficients[entry] != value for entry, value in forced.items()):
        raise AssertionError("coarse map overwrote a forced physical orbit")

    target = coarse_open_completion_target()
    if not np.isclose(target.coefficient, 1 / 32, atol=0):
        raise AssertionError(("coarse target coefficient", target))
    if not np.isclose(target.optimal_beta, 0.779698447178, atol=3e-9):
        raise AssertionError(("coarse target beta", target))
    if not np.isclose(target.optimized_total, 0.322669154028, atol=3e-11):
        raise AssertionError(("coarse target total", target))
    if target.threshold_slack <= 0.0106:
        raise AssertionError(("coarse target slack", target))

    print(
        "repaired open-profile budget passed: "
        f"open_profiles={len(profiles)},"
        f"repaired_baseline={repaired.baseline_total:.12g},"
        f"extended_common={extended.common_coefficient_threshold:.12g},"
        f"forced_entries={target.forced_entries},"
        f"coarse_entries={target.coarse_entries},"
        f"coarse_beta={target.optimal_beta:.12g},"
        f"coarse_total={target.optimized_total:.12g},"
        f"coarse_slack={target.threshold_slack:.12g}"
    )


if __name__ == "__main__":
    main()
