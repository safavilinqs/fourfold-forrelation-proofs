#!/usr/bin/env python3
"""Reproduce the quantitative q64 active-platform feasibility screen.

The mathematical active threshold is exact.  The 98 percent detector value
is a published benchmark, not an assumption in the separation theorem and
not a certification of a complete device.  The remaining calculations show
what optical budget would remain if that benchmark were attained.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, getcontext
from json import dumps
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
getcontext().prec = 80

ONE = Decimal(1)
TEN = Decimal(10)
DETECTOR_BENCHMARK = Decimal("0.98")
DIMENSION = 4096
SIGN_MODES = 4 * DIMENSION
BUTTERFLY_STAGES_PER_HADAMARD = 12
GEOMETRIC_MEAN_BUTTERFLY_DEPTH = 18


def threshold_expectation() -> Decimal:
    """Return the positive root of 3 mu^3 - 9 mu + 2 = 0 below 1/4."""

    value = Decimal("0.226")
    for _ in range(30):
        residual = 3 * value**3 - 9 * value + 2
        derivative = 9 * value**2 - 9
        value -= residual / derivative
    if not (Decimal("0.22") < value < Decimal("0.25")):
        raise AssertionError(value)
    return value


def db_loss(transmission: Decimal) -> Decimal:
    """Return minus ten log10 of a power transmission."""

    if not 0 < transmission <= 1:
        raise ValueError(transmission)
    return -TEN * transmission.ln() / TEN.ln()


def diagnostic() -> dict[str, object]:
    mu = threshold_expectation()
    minimum_contrast = 4 * mu
    residual_contrast = minimum_contrast / DETECTOR_BENCHMARK
    per_stage_transmission = (
        residual_contrast.ln() / GEOMETRIC_MEAN_BUTTERFLY_DEPTH
    ).exp()
    return {
        "schema": "round4_q64_experimental_feasibility_gate_v1",
        "verdict": "NOT_YET_EXPERIMENTALLY_CREDIBLE",
        "scope": (
            "engineering screen for the active implementation; not part of "
            "the passive lower-bound theorem and not an impossibility theorem"
        ),
        "mathematical_requirements": {
            "dimension": DIMENSION,
            "sign_modes": SIGN_MODES,
            "minimum_flag_expectation": str(mu),
            "minimum_combined_contrast": str(minimum_contrast),
            "butterfly_stages_per_hadamard": BUTTERFLY_STAGES_PER_HADAMARD,
            "left_hadamards": 1,
            "right_hadamards": 2,
        },
        "nominal_detector_benchmark_screen": {
            "detector_system_efficiency": str(DETECTOR_BENCHMARK),
            "minimum_remaining_nondetector_contrast": str(residual_contrast),
            "maximum_geometric_mean_nondetector_loss_db": str(
                db_loss(residual_contrast)
            ),
            "geometric_mean_butterfly_depth": GEOMETRIC_MEAN_BUTTERFLY_DEPTH,
            "minimum_identical_stage_power_transmission": str(
                per_stage_transmission
            ),
            "maximum_identical_stage_loss_db": str(
                db_loss(per_stage_transmission)
            ),
            "interpretation": (
                "illustrative best-case allocation with perfect masks, "
                "routing, overlap, phase, and zero additive bias"
            ),
        },
        "primary_evidence_snapshot": [
            {
                "platform": "SNSPD",
                "demonstration": "98.0 +/- 0.5 percent system efficiency at 1550 nm",
                "url": "https://doi.org/10.1364/OPTICA.400751",
            },
            {
                "platform": "ultrafast time bins",
                "demonstration": (
                    "362 programmable unitaries through dimension 8 and a "
                    "passive network through 36 modes, with fidelity above 97 percent"
                ),
                "url": "https://doi.org/10.1103/PhysRevLett.133.090601",
            },
            {
                "platform": "time-domain Gaussian photonics",
                "demonstration": "216-mode nonuniversal three-loop processor",
                "url": "https://doi.org/10.1038/s41586-022-04725-x",
            },
            {
                "platform": "complex-medium spatial optics",
                "demonstration": (
                    "approximately 200 mixer modes; characterized Fourier "
                    "circuits only through dimension 7"
                ),
                "url": "https://doi.org/10.1038/s41567-023-02319-6",
            },
        ],
    }


def artifact_text(result: dict[str, object]) -> str:
    return dumps(result, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-artifact", action="store_true")
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_experimental_feasibility_gate.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    screen = result["nominal_detector_benchmark_screen"]
    assert isinstance(screen, dict)
    print(
        "q64 experimental feasibility gate: "
        f"contrast={result['mathematical_requirements']['minimum_combined_contrast']},"
        f"residual_loss_db={screen['maximum_geometric_mean_nondetector_loss_db']},"
        f"stage_transmission={screen['minimum_identical_stage_power_transmission']},"
        "verdict=not_yet_experimentally_credible"
    )


if __name__ == "__main__":
    main()
