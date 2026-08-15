#!/usr/bin/env python3
"""Cross-check the concise q64 theorem and experimental closeout package."""

from __future__ import annotations

from decimal import Decimal, getcontext
from json import loads
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
getcontext().prec = 90


def load_artifact(name: str) -> dict[str, object]:
    return loads((ROOT / "artifacts" / name).read_text(encoding="utf-8"))


def main() -> None:
    adaptive = load_artifact("q64_adaptive_tree_frontier.json")["result"]
    active = load_artifact("active_six_resource_row.json")["result"]
    feasibility = load_artifact("q64_experimental_feasibility_gate.json")

    if (adaptive["dimension"], adaptive["sign_modes"]) != (4096, 16384):
        raise AssertionError(("adaptive size", adaptive))
    if (active["dimension"], active["sign_modes"]) != (4096, 16384):
        raise AssertionError(("active size", active))
    if (active["total_hard_dose"], active["majority_error_exact"]) != (6, "81/256"):
        raise AssertionError(("active theorem", active))
    if adaptive["frontier_multiplier_numerator"] != adaptive[
        "frontier_multiplier_denominator"
    ]:
        raise AssertionError(("adaptive multiplier", adaptive))

    total = Decimal(adaptive["adaptive_total_upper"])
    threshold = Decimal(adaptive["threshold_lower"])
    reserve = Decimal(adaptive["reserve_margin_lower"])
    if not reserve <= threshold - total < reserve + Decimal("1e-75"):
        raise AssertionError(("adaptive reserve", threshold - total, reserve))
    if (Decimal(1) - total) / 2 <= Decimal(1) / 3:
        raise AssertionError(("passive average error", total))
    if feasibility["verdict"] != "NOT_YET_EXPERIMENTALLY_CREDIBLE":
        raise AssertionError(("feasibility verdict", feasibility))

    theorem_path = ROOT / "notes" / "PAPER_CLOSEOUT_THEOREM_PACKAGE.md"
    theorem = theorem_path.read_text(encoding="utf-8")
    for phrase in (
        "STATUS: NUMBER-SECTOR-INCOHERENT MATHEMATICAL SEPARATION CERTIFIED AT $N=4096$",
        "EXPERIMENTAL LABEL: IDEAL RESOURCE SPECIFICATION",
        "F_{4,H}(x)\\ge\\frac14",
        "\\beta=19/25",
        "0.2609692247922079249341809573938165614",
        "D_{\\mathsf P}^{\\rm hard}(4096)>6",
        "D_{\\mathsf A}^{\\rm hard}(4096)\\le6",
        "81}{256",
        "block diagonal in total signal photon number",
        "extension to fresh probes with coherence between different total signal-number sectors is open",
        "272 unbalanced high-sector split/state incidences",
        "./run_round4_checks.sh",
    ):
        if phrase not in theorem:
            raise AssertionError(("paper theorem contract", phrase))

    for relative in (
        "notes/PAPER_CLOSEOUT_THEOREM_PACKAGE.md",
        "notes/EXPERIMENTAL_FEASIBILITY_DECISION.md",
    ):
        text = (ROOT / relative).read_text(encoding="utf-8")
        if "\\(" in text or "\\[" in text:
            raise AssertionError(("non-Obsidian delimiter", relative))

    print(
        "q64 paper closeout package passed: "
        f"adaptive_total={total},active_error=81/256,"
        "hardware=not_yet_experimentally_credible"
    )


if __name__ == "__main__":
    main()
