#!/usr/bin/env python3
"""Regression for the q64 active-platform feasibility arithmetic."""

from __future__ import annotations

from decimal import Decimal
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_experimental_feasibility_gate import (  # noqa: E402
    artifact_text,
    diagnostic,
)


def main() -> None:
    result = diagnostic()
    requirements = result["mathematical_requirements"]
    screen = result["nominal_detector_benchmark_screen"]
    if result["verdict"] != "NOT_YET_EXPERIMENTALLY_CREDIBLE":
        raise AssertionError(result["verdict"])
    if (requirements["dimension"], requirements["sign_modes"]) != (4096, 16384):
        raise AssertionError(requirements)

    contrast = Decimal(requirements["minimum_combined_contrast"])
    residual = Decimal(screen["minimum_remaining_nondetector_contrast"])
    stage = Decimal(screen["minimum_identical_stage_power_transmission"])
    if abs(contrast - Decimal("0.904294855156833")) > Decimal("3e-15"):
        raise AssertionError(("contrast", contrast))
    if residual != contrast / Decimal("0.98"):
        raise AssertionError(("residual", residual))
    if abs(stage**18 - residual) > Decimal("1e-70"):
        raise AssertionError(("stage depth", stage**18, residual))
    if not Decimal("0.34") < Decimal(
        screen["maximum_geometric_mean_nondetector_loss_db"]
    ) < Decimal("0.36"):
        raise AssertionError(("loss budget", screen))

    committed = (
        ROOT / "artifacts" / "q64_experimental_feasibility_gate.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 experimental feasibility artifact")

    decision = (ROOT / "notes" / "EXPERIMENTAL_FEASIBILITY_DECISION.md").read_text(
        encoding="utf-8"
    )
    for phrase in (
        "VERDICT: NOT YET EXPERIMENTALLY CREDIBLE",
        "$4096$-dimensional coherent Sylvester transform",
        "0.35 dB",
        "not a physical impossibility theorem",
    ):
        if phrase not in decision:
            raise AssertionError(("experimental decision contract", phrase))

    print(
        "q64 experimental feasibility gate passed: "
        f"residual_loss_db={screen['maximum_geometric_mean_nondetector_loss_db']},"
        f"stage_transmission={stage},"
        "verdict=not_yet_experimentally_credible"
    )


if __name__ == "__main__":
    main()
