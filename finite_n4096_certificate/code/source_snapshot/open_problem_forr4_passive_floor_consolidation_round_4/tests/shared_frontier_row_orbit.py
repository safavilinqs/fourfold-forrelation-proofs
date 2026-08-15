#!/usr/bin/env python3
"""Regression for the native-q32 shared-frontier row-orbit screen."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from shared_frontier_row_orbit import artifact_text, diagnostic  # noqa: E402


def main() -> None:
    result = diagnostic()
    observed_discrete = (
        result.order,
        result.basis_dimension,
        result.occupation_states,
        result.translation_orbit_size,
        result.searched_trials,
        result.selected_trial,
        result.relevant_occupation_edges,
        result.activated_orbits,
        result.activated_entries,
        result.matching_pairs,
        result.nonzero_pairs,
    )
    expected_discrete = (32, 960, 30, 32, 9, 8, 67, 6, 6, 50_112, 1_984)
    if observed_discrete != expected_discrete:
        raise AssertionError(("row-orbit discrete result", observed_discrete))

    observed = (
        result.attenuated_separate_nuclear,
        result.attenuated_joint_nuclear,
        result.attenuated_cancellation_ratio,
        result.current_frontier_perron_contribution,
        result.separate_to_current_frontier_ratio,
    )
    expected = (
        4.794169961088444e-13,
        4.518664344625811e-13,
        0.9425331978843562,
        0.028226470553668472,
        1.6984659672462544e-11,
    )
    if not np.allclose(observed, expected, rtol=3e-11, atol=2e-25):
        raise AssertionError(("row-orbit numeric result", observed))
    if result.activated_orbits >= 10:
        raise AssertionError("row-orbit parity screen unexpectedly became broad")
    if result.separate_to_current_frontier_ratio >= 1e-9:
        raise AssertionError("row-orbit screen unexpectedly became competitive")

    committed = (ROOT / "artifacts" / "shared_frontier_row_orbit.json").read_text(
        encoding="utf-8"
    )
    if committed != artifact_text(result):
        raise AssertionError("stale shared-frontier row-orbit artifact")

    print(
        "shared frontier row orbit passed: "
        f"basis={result.basis_dimension},"
        f"orbits={result.activated_orbits},"
        f"nonzero_pairs={result.nonzero_pairs},"
        f"separate={result.attenuated_separate_nuclear:.12g},"
        f"joint={result.attenuated_joint_nuclear:.12g},"
        f"frontier_ratio={result.separate_to_current_frontier_ratio:.12g},"
        "decision=close_row_translation_lower_witness"
    )


if __name__ == "__main__":
    main()
