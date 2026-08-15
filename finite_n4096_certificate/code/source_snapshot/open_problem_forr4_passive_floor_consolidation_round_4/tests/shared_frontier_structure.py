#!/usr/bin/env python3
"""Regression for the exact leading-frontier occupation reduction."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from shared_frontier_structure import artifact_text, audit  # noqa: E402


def main() -> None:
    result = audit()
    frontier = result["frontier"]
    kernel = result["kernel"]
    layers = result["degree_layers"]
    patterns = result["profile_patterns"]

    if not isinstance(frontier, dict) or not isinstance(kernel, dict):
        raise AssertionError("malformed frontier structure")
    observed_frontier = (
        frontier["orbits"],
        frontier["entries"],
        frontier["unordered_profile_patterns"],
        kernel["occupation_states"],
        kernel["occupation_edges"],
        kernel["connected_components"],
        kernel["cross_degree_shared_edges"],
    )
    if observed_frontier != (51, 198, 6, 125, 241, 8, 25):
        raise AssertionError(("frontier reduction", observed_frontier))
    if not np.isclose(
        frontier["fraction_of_unresolved_perron_contribution"],
        0.9054673281988463,
        rtol=0,
        atol=2e-15,
    ):
        raise AssertionError(("frontier fraction", frontier))

    expected_layers = (
        (10, 5, 100, 500, {"0": 100, "1": 400}, [5, 6], 217),
        (12, 6, 98, 98, {"0": 98}, [6], 49),
    )
    observed_layers = tuple(
        (
            row["degree"],
            row["balanced_split_degree"],
            row["entries"],
            row["compatible_terms"],
            row["shared_intersection_term_counts"],
            row["occupation_layers"],
            row["unique_occupation_edges"],
        )
        for row in layers
    )
    if observed_layers != expected_layers:
        raise AssertionError(("degree layers", observed_layers))

    expected_components = (
        (5, 11),
        (5, 11),
        (5, 11),
        (5, 11),
        (6, 24),
        (6, 20),
        (6, 19),
        (6, 18),
    )
    observed_components = tuple(
        (row["occupation_layer"], row["states"]) for row in kernel["components"]
    )
    if observed_components != expected_components:
        raise AssertionError(("components", observed_components))

    expected_patterns = (
        ((5, 3, 1, 1), 21, 84, 189),
        ((5, 3, 3, 1), 15, 60, 30),
        ((7, 3, 1, 1), 4, 16, 8),
        ((7, 1, 1, 1), 4, 16, 40),
        ((5, 5, 1, 1), 5, 16, 8),
        ((3, 3, 3, 3), 2, 6, 3),
    )
    observed_patterns = tuple(
        (
            tuple(row["sorted_profile"]),
            row["orbits"],
            row["entries"],
            row["occupation_edges"],
        )
        for row in patterns
    )
    if observed_patterns != expected_patterns:
        raise AssertionError(("patterns", observed_patterns))

    committed = (ROOT / "artifacts" / "shared_frontier_structure.json").read_text(
        encoding="utf-8"
    )
    if committed != artifact_text(result):
        raise AssertionError("stale shared-frontier structure artifact")

    print(
        "shared frontier structure passed: "
        "orbits=51,entries=198,patterns=6,states=125,edges=241,components=8"
    )


if __name__ == "__main__":
    main()
