#!/usr/bin/env python3
"""Independent checks for the q64 adaptive-lift arithmetic and interface claims."""

from __future__ import annotations

from fractions import Fraction
from json import loads
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_adaptive_interface_gate import artifact_text, diagnostic  # noqa: E402


def direct_posterior_counterexample() -> None:
    """Reject reuse of an unconditional child bound after conditioning."""

    plus = {(-1, -1): Fraction(1, 2), (1, 1): Fraction(1, 2)}
    minus = {(-1, 1): Fraction(1, 2), (1, -1): Fraction(1, 2)}

    def marginal(distribution: dict[tuple[int, int], Fraction], axis: int) -> dict[int, Fraction]:
        result = {-1: Fraction(0), 1: Fraction(0)}
        for point, mass in distribution.items():
            result[point[axis]] += mass
        return result

    root_plus = marginal(plus, 0)
    root_minus = marginal(minus, 0)
    root_tv = sum(abs(root_plus[x] - root_minus[x]) for x in (-1, 1)) / 2
    if root_tv != 0:
        raise AssertionError("toy root should carry no hypothesis information")

    for observed in (-1, 1):
        child_plus = {
            point[1]: mass / root_plus[observed]
            for point, mass in plus.items()
            if point[0] == observed
        }
        child_minus = {
            point[1]: mass / root_minus[observed]
            for point, mass in minus.items()
            if point[0] == observed
        }
        child_tv = sum(
            abs(child_plus.get(x, 0) - child_minus.get(x, 0))
            for x in (-1, 1)
        ) / 2
        if child_tv != 1:
            raise AssertionError((observed, child_tv))

    joint_tv = sum(
        abs(plus.get(point, 0) - minus.get(point, 0))
        for point in set(plus) | set(minus)
    ) / 2
    if joint_tv != 1:
        raise AssertionError(joint_tv)


def schur_amplification_checks() -> None:
    """Stress the factorization bound after adding auxiliary outcome vectors."""

    generator = np.random.default_rng(20260718)
    for _ in range(80):
        rows, columns, rank, auxiliary = 5, 6, 4, 3
        left = generator.normal(size=(rows, rank))
        right = generator.normal(size=(columns, rank))
        left /= np.maximum(np.linalg.norm(left, axis=1, keepdims=True), 1.0)
        right /= np.maximum(np.linalg.norm(right, axis=1, keepdims=True), 1.0)
        kernel = left @ right.T

        outcome_left = generator.normal(size=(rows, auxiliary))
        outcome_right = generator.normal(size=(columns, auxiliary))
        outcome_left /= np.maximum(
            np.linalg.norm(outcome_left, axis=1, keepdims=True), 1.0
        )
        outcome_right /= np.maximum(
            np.linalg.norm(outcome_right, axis=1, keepdims=True), 1.0
        )
        amplified = kernel * (outcome_left @ outcome_right.T)

        p = generator.dirichlet(np.ones(rows))
        w = generator.dirichlet(np.ones(columns))
        weighted = np.sqrt(p)[:, None] * amplified * np.sqrt(w)[None, :]
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > 1 + 2e-12:
            raise AssertionError(("Schur amplification exceeded factorization", nuclear))


def main() -> None:
    result = diagnostic()
    additive = Fraction(result.additive_overhead_cap)
    occurrence_cap = Fraction(result.occurrence_multiplier_cap)
    whole_cap = Fraction(result.whole_total_multiplier_cap)
    if not (Fraction(6347, 100000) < additive < Fraction(6348, 100000)):
        raise AssertionError(("additive cap", additive))
    if not (Fraction(1238, 1000) < occurrence_cap < Fraction(1239, 1000)):
        raise AssertionError(("occurrence cap", occurrence_cap))
    if not (Fraction(1236, 1000) < whole_cap < Fraction(1237, 1000)):
        raise AssertionError(("whole cap", whole_cap))
    if not result.zero_overhead_passes or not result.six_fifths_occurrence_multiplier_passes:
        raise AssertionError("known passing multipliers")
    if result.five_fourths_occurrence_multiplier_passes or result.twofold_occurrence_multiplier_passes:
        raise AssertionError("known failing multipliers")
    if result.verdict != "INCOMPLETE":
        raise AssertionError(result.verdict)

    direct_posterior_counterexample()
    schur_amplification_checks()

    committed = (ROOT / "artifacts" / "q64_adaptive_interface_gate.json").read_text(
        encoding="utf-8"
    )
    if committed != artifact_text(result):
        raise AssertionError("stale adaptive-interface artifact")
    payload = loads(committed)
    if payload["result"]["verdict"] != "INCOMPLETE":
        raise AssertionError("artifact verdict")
    print(
        "q64 adaptive interface gate passed: "
        f"additive_cap={result.additive_overhead_cap},"
        f"occurrence_multiplier_cap={result.occurrence_multiplier_cap},"
        "posterior_reuse=rejected,schur_amplification=passed,"
        "verdict=incomplete"
    )


if __name__ == "__main__":
    main()
