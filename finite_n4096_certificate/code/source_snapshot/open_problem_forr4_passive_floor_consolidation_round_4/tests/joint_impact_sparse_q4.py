#!/usr/bin/env python3
"""Regression for the compatible leading-frontier physical law."""

from __future__ import annotations

from math import comb, prod
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from joint_impact_sparse_q4 import (  # noqa: E402
    FRONTIER_ORBITS,
    artifact_text,
    diagnostics,
    frontier,
    occupation_law,
)


def main() -> None:
    orbit_entries = frontier()
    if len(orbit_entries) != FRONTIER_ORBITS:
        raise AssertionError(("frontier", len(orbit_entries)))
    if any(max(profile) <= 4 for orbit in orbit_entries for profile, _ in orbit):
        raise AssertionError("a leading orbit was unexpectedly q=2-realizable")

    _, law = occupation_law()
    full_q4_dimension = sum(
        prod(comb(16, occupation) for occupation in state) for state, _ in law
    )
    if full_q4_dimension != 18_904_064:
        raise AssertionError(("full q4 basis", full_q4_dimension))

    q4, q32 = diagnostics()
    expected_discrete = (
        q4.basis_dimension,
        q4.support_states,
        q4.selected_trial,
        q4.activated_orbits,
        q4.activated_entries,
        q4.nonzero_matrix_entries,
        q32.activated_orbits,
        q32.activated_entries,
        q32.nonzero_matrix_entries,
    )
    if expected_discrete != (180, 30, 16, 16, 26, 128, 16, 26, 128):
        raise AssertionError(("joint discrete result", expected_discrete))

    observed = (
        q4.separate_nuclear,
        q4.joint_nuclear,
        q4.attenuated_separate_nuclear,
        q4.attenuated_joint_nuclear,
        q32.separate_nuclear,
        q32.joint_nuclear,
        q32.attenuated_separate_nuclear,
        q32.attenuated_joint_nuclear,
    )
    expected = (
        0.0026244357009079368,
        0.0018718547844748336,
        0.00019583215915551996,
        0.0001461225757568427,
        1.3658615751792888e-08,
        1.1571372011432606e-08,
        1.1152267139576358e-09,
        9.52512381241109e-10,
    )
    if not np.allclose(observed, expected, rtol=2e-11, atol=2e-16):
        raise AssertionError(("joint norm result", observed))
    if not 0 < q4.attenuated_cancellation_ratio < 0.8:
        raise AssertionError(("q4 cancellation", q4))
    if not 0 < q32.attenuated_cancellation_ratio < 0.9:
        raise AssertionError(("q32 cancellation", q32))

    committed = (ROOT / "artifacts" / "joint_impact_sparse_diagnostic.json").read_text(
        encoding="utf-8"
    )
    regenerated = artifact_text(q4, q32)
    if committed != regenerated:
        raise AssertionError("stale joint physical-law artifact")

    print(
        "joint impact sparse law passed: "
        f"q2_inaccessible_orbits={FRONTIER_ORBITS},"
        f"full_q4_basis={full_q4_dimension},"
        f"sampled_basis={q4.basis_dimension},"
        f"q4_orbits={q4.activated_orbits},"
        f"q4_ratio={q4.attenuated_cancellation_ratio:.12g},"
        f"q32_orbits={q32.activated_orbits},"
        f"q32_ratio={q32.attenuated_cancellation_ratio:.12g}"
    )


if __name__ == "__main__":
    main()
