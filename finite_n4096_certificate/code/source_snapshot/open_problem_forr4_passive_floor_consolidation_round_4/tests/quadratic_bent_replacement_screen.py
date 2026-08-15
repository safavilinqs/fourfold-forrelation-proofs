#!/usr/bin/env python3
"""Regression for the quadratic-bent finite-size replacement screen."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from quadratic_bent_replacement_screen import (  # noqa: E402
    artifact_text,
    screen,
)


def main() -> None:
    result = screen()
    counts = (
        result.nonminimal_profiles,
        result.raw_profile_splits,
        result.compatible_profile_splits,
        result.compatible_by_degree,
    )
    expected_counts = (69, 7904, 2284, {6: 128, 8: 476, 10: 920, 12: 760})
    if counts != expected_counts:
        raise AssertionError(("replacement inventory", counts))

    observed = (
        result.optimistic_zero_higher_floor,
        result.optimistic_floor_margin,
        result.common_coefficient_gate,
        result.common_gate_times_dimension,
        result.half_over_dimension_total,
        result.half_over_dimension_margin,
        result.one_over_dimension_total,
        result.one_over_dimension_overshoot,
        result.degree_only_gates[6],
        result.degree_only_gates[8],
        result.degree_only_gates[10],
        result.degree_only_gates[12],
        result.known_endpoint_five_one_coefficient,
    )
    expected = (
        0.28151203289098925,
        0.051821300442344065,
        0.0005298085951150572,
        0.5425240013978185,
        0.3286702351273519,
        0.004663098205981431,
        0.3876188908534863,
        0.05428555752015296,
        0.0016204016251463356,
        0.0013280986879687956,
        0.002295000996958223,
        0.009676567754345946,
        1 / 511,
    )
    if not np.allclose(observed, expected, rtol=3e-10, atol=3e-12):
        raise AssertionError(("replacement gates", observed))
    if result.half_over_dimension_total >= 1 / 3:
        raise AssertionError(("half-over-N target", result))
    if result.one_over_dimension_total <= 1 / 3:
        raise AssertionError(("one-over-N rejection", result))
    if result.decision != "not_promoted_current_signed_permutation_remains_lead":
        raise AssertionError(("replacement decision", result.decision))

    committed = (
        ROOT / "artifacts" / "quadratic_bent_replacement_screen.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale quadratic-bent replacement artifact")

    print(
        "quadratic-bent replacement screen passed: "
        f"compatible={result.compatible_profile_splits},"
        f"floor={result.optimistic_zero_higher_floor:.12g},"
        f"common_gate={result.common_coefficient_gate:.12g},"
        f"half_over_N_margin={result.half_over_dimension_margin:.12g},"
        f"one_over_N_overshoot={result.one_over_dimension_overshoot:.12g},"
        f"decision={result.decision}"
    )


if __name__ == "__main__":
    main()
