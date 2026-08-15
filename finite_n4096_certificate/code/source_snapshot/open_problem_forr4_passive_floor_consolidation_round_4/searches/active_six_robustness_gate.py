#!/usr/bin/env python3
"""Exact active-side contrast gate for the three-flag dose-six protocol.

Assign an independent fair sign to every no-click or erased flag.  If the
signed expectation of each retained binary flag at the promise boundary is
``mu``, majority decoding three independent flags has error

    ((1 - mu) / 2)^2 * (2 + mu).

This module solves the exact one-third gate and translates it into a single
coherent-throughput/visibility requirement for the active resource row.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from math import acos, sqrt
from pathlib import Path

from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
PROMISE_MAGNITUDE = 1 / 4
ERROR_THRESHOLD = 1 / 3
IDEAL_FLAG_EXPECTATION = PROMISE_MAGNITUDE


@dataclass(frozen=True)
class ActiveRobustnessGate:
    flags: int
    hard_dose: int
    ideal_flag_expectation: float
    ideal_majority_error: float
    threshold_flag_expectation: float
    threshold_single_flag_correct_probability: float
    minimum_multiplicative_contrast: float
    maximum_additive_expectation_error_at_unit_contrast: float
    equal_two_pass_minimum_power_transmission_per_pass: float
    phase_only_maximum_absolute_radians: float
    phase_only_maximum_absolute_degrees: float


def majority_error(expectation: float) -> float:
    if not -1 <= expectation <= 1:
        raise ValueError(("binary expectation outside [-1,1]", expectation))
    wrong = (1 - expectation) / 2
    return wrong * wrong * (2 + expectation)


def diagnostic() -> ActiveRobustnessGate:
    threshold_expectation = float(
        brentq(
            lambda value: majority_error(value) - ERROR_THRESHOLD,
            0.0,
            IDEAL_FLAG_EXPECTATION,
            xtol=1e-15,
        )
    )
    contrast = threshold_expectation / PROMISE_MAGNITUDE
    phase = acos(contrast)
    return ActiveRobustnessGate(
        flags=3,
        hard_dose=6,
        ideal_flag_expectation=IDEAL_FLAG_EXPECTATION,
        ideal_majority_error=majority_error(IDEAL_FLAG_EXPECTATION),
        threshold_flag_expectation=threshold_expectation,
        threshold_single_flag_correct_probability=(1 + threshold_expectation) / 2,
        minimum_multiplicative_contrast=contrast,
        maximum_additive_expectation_error_at_unit_contrast=(
            IDEAL_FLAG_EXPECTATION - threshold_expectation
        ),
        equal_two_pass_minimum_power_transmission_per_pass=sqrt(contrast),
        phase_only_maximum_absolute_radians=phase,
        phase_only_maximum_absolute_degrees=phase * 180 / 3.141592653589793,
    )


def artifact_text(result: ActiveRobustnessGate) -> str:
    payload = {
        "schema": "round4_active_six_robustness_gate_v1",
        "result": asdict(result),
        "physical_interpretation": (
            "if g is the product of detection efficiency, geometric-mean "
            "arm power transmission, mode overlap, and differential-phase "
            "cosine, and b is worst signed additive flag bias, require "
            "g/4 - b above the threshold flag expectation"
        ),
        "evidence_label": (
            "exact three-independent-flag majority calculation; scalar "
            "imperfection translation; not a device-specific noise model"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic active robustness gate",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "active_six_robustness_gate.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "active six-dose robustness gate: "
        f"threshold_expectation={result.threshold_flag_expectation:.12g},"
        f"minimum_contrast={result.minimum_multiplicative_contrast:.12g},"
        "equal_two_pass_minimum="
        f"{result.equal_two_pass_minimum_power_transmission_per_pass:.12g},"
        f"phase_only_max_degrees={result.phase_only_maximum_absolute_degrees:.12g},"
        "status=exact_scalar_active_gate"
    )


if __name__ == "__main__":
    main()
