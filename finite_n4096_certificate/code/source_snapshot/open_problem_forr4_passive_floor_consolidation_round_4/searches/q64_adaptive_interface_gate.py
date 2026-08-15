#!/usr/bin/env python3
"""Exact arithmetic gate for the q64 adaptive-lift interface."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from json import dumps, loads
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "artifacts" / "q64_complete_outward_ledger.json"


@dataclass(frozen=True)
class AdaptiveInterfaceGate:
    perron_upper: str
    promise_upper: str
    one_batch_total_upper: str
    threshold_lower: str
    additive_overhead_cap: str
    occurrence_multiplier_cap: str
    whole_total_multiplier_cap: str
    zero_overhead_passes: bool
    six_fifths_occurrence_multiplier_passes: bool
    five_fourths_occurrence_multiplier_passes: bool
    twofold_occurrence_multiplier_passes: bool
    verdict: str
    missing_lemma: str


def _fraction(text: str) -> Fraction:
    return Fraction(text)


def _decimal_floor(value: Fraction, digits: int = 60) -> str:
    scale = 10**digits
    units = value.numerator * scale // value.denominator
    integer, remainder = divmod(units, scale)
    return f"{integer}.{remainder:0{digits}d}"


def _load_ledger() -> dict[str, object]:
    payload = loads(LEDGER.read_text(encoding="utf-8"))
    if payload.get("schema") != "round4_q64_complete_outward_ledger_v3":
        raise AssertionError("unexpected complete-ledger schema")
    return payload["result"]


def diagnostic() -> AdaptiveInterfaceGate:
    result = _load_ledger()
    perron = _fraction(str(result["collatz_perron_upper"]))
    promise = _fraction(str(result["promise_loss_upper"]))
    total = _fraction(str(result["total_upper"]))
    threshold = _fraction(str(result["reserve_threshold"]))
    if perron + promise > total:
        raise AssertionError("complete total does not dominate its components")

    additive = threshold - total
    occurrence_cap = (threshold - promise) / perron
    whole_cap = threshold / total

    def occurrence_passes(multiplier: Fraction) -> bool:
        return multiplier * perron + promise < threshold

    return AdaptiveInterfaceGate(
        perron_upper=str(result["collatz_perron_upper"]),
        promise_upper=str(result["promise_loss_upper"]),
        one_batch_total_upper=str(result["total_upper"]),
        threshold_lower=str(result["reserve_threshold"]),
        additive_overhead_cap=_decimal_floor(additive),
        occurrence_multiplier_cap=_decimal_floor(occurrence_cap),
        whole_total_multiplier_cap=_decimal_floor(whole_cap),
        zero_overhead_passes=total < threshold,
        six_fifths_occurrence_multiplier_passes=occurrence_passes(Fraction(6, 5)),
        five_fourths_occurrence_multiplier_passes=occurrence_passes(Fraction(5, 4)),
        twofold_occurrence_multiplier_passes=occurrence_passes(Fraction(2)),
        verdict="INCOMPLETE",
        missing_lemma=(
            "One globally normalized operator-valued occurrence frontier for "
            "every posterior-selected passive tree under branchwise hard dose six."
        ),
    )


def artifact_text(result: AdaptiveInterfaceGate) -> str:
    payload = {
        "schema": "round4_q64_adaptive_interface_gate_v1",
        "result": asdict(result),
        "interpretation": (
            "The occurrence multiplier applies only to the certified Perron term; "
            "the promise loss is retained once. A candidate multiplier must be "
            "strictly below the displayed inward-rounded cap."
        ),
        "scope": (
            "Arithmetic acceptance gate and interface audit only. This artifact "
            "does not prove the adaptive frontier-normalization lemma."
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.output is not None:
        arguments.output.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 adaptive interface gate: "
        f"additive_cap={result.additive_overhead_cap},"
        f"occurrence_multiplier_cap={result.occurrence_multiplier_cap},"
        f"whole_total_multiplier_cap={result.whole_total_multiplier_cap},"
        f"verdict={result.verdict.lower()}"
    )


if __name__ == "__main__":
    main()
