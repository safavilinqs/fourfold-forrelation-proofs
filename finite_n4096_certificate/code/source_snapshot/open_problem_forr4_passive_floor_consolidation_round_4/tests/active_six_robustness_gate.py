#!/usr/bin/env python3
"""Regression for the exact scalar active robustness gate."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from active_six_robustness_gate import (  # noqa: E402
    artifact_text,
    diagnostic,
    majority_error,
)


def main() -> None:
    result = diagnostic()
    if (result.flags, result.hard_dose) != (3, 6):
        raise AssertionError(("active robustness resources", result))
    observed = (
        result.ideal_flag_expectation,
        result.ideal_majority_error,
        result.threshold_flag_expectation,
        result.threshold_single_flag_correct_probability,
        result.minimum_multiplicative_contrast,
        result.maximum_additive_expectation_error_at_unit_contrast,
        result.equal_two_pass_minimum_power_transmission_per_pass,
        result.phase_only_maximum_absolute_radians,
        result.phase_only_maximum_absolute_degrees,
    )
    expected = (
        0.25,
        81 / 256,
        0.22607371378920826,
        0.6130368568946041,
        0.904294855156833,
        0.023926286210791736,
        0.950944191399702,
        0.4410712531917783,
        25.27152127243504,
    )
    if not np.allclose(observed, expected, rtol=3e-12, atol=3e-14):
        raise AssertionError(("active robustness numeric result", observed))
    if not np.isclose(
        majority_error(result.threshold_flag_expectation), 1 / 3, atol=1e-14
    ):
        raise AssertionError("active robustness threshold does not reach 1/3")

    committed = (ROOT / "artifacts" / "active_six_robustness_gate.json").read_text(
        encoding="utf-8"
    )
    if committed != artifact_text(result):
        raise AssertionError("stale active robustness artifact")

    print(
        "active six-dose robustness gate passed: "
        f"threshold_expectation={result.threshold_flag_expectation:.12g},"
        f"minimum_contrast={result.minimum_multiplicative_contrast:.12g},"
        "equal_two_pass_minimum="
        f"{result.equal_two_pass_minimum_power_transmission_per_pass:.12g},"
        "status=exact_scalar_active_gate"
    )


if __name__ == "__main__":
    main()
