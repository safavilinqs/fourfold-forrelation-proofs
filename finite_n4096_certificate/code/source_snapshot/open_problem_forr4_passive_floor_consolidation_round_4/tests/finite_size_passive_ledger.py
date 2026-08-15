#!/usr/bin/env python3
"""Regression for the generated N=1024 finite-size routing ledger."""

from __future__ import annotations

from json import loads
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from finite_size_passive_ledger import build_audit, write_artifacts  # noqa: E402


def main() -> None:
    audit = build_audit()
    expected_counts = {
        "open_profile_splits": 6016,
        "off_face_generic_certified_at_1_over_q": 1138,
        "balanced_entries": 888,
        "proved_arbitrary_upper_entries": 40,
        "physical_lower_witness_substitutions": 8,
        "local_slice_substitutions": 4,
        "provisional_common_target_entries": 836,
        "unresolved_entries": 848,
        "unresolved_orbits": 224,
    }
    if audit.counts != expected_counts:
        raise AssertionError(("ledger counts", audit.counts))

    expected_frontier = (
        0.13626448896970556,
        0.510656928241998,
        0.905467328198846,
        16,
        51,
    )
    observed_frontier = (
        audit.impact_frontier["leading_orbit_fraction"],
        audit.impact_frontier["leading_16_orbits_fraction"],
        audit.impact_frontier["leading_51_orbits_fraction"],
        audit.impact_frontier["orbits_reaching_50_percent"],
        audit.impact_frontier["orbits_reaching_90_percent"],
    )
    if not np.allclose(observed_frontier, expected_frontier, atol=5e-14):
        raise AssertionError(("impact frontier", observed_frontier))

    expected_totals = (
        0.248224749551,
        0.279758546919,
        0.322669154028,
        0.779315347679,
        0.333132605485488,
        0.000200727847845,
    )
    observed_totals = (
        audit.totals.degree_eight_total,
        audit.totals.known_high_total,
        audit.totals.coarse_before_ten_total,
        audit.totals.ten_theorem_beta,
        audit.totals.ten_theorem_total,
        audit.totals.ten_theorem_margin,
    )
    if not np.allclose(observed_totals, expected_totals, atol=4e-10):
        raise AssertionError(("ledger totals", observed_totals))

    coefficients = tuple(family.coefficient for family in audit.theorem_families)
    expected_coefficients = (
        0.0934752745775,
        0.0162724692796,
        0.173742800847,
        0.0250967461185,
        0.0311889051224,
        0.0422410016249,
        0.0370952793157,
        0.0285281522923,
        0.0250919471547,
        0.0462425962446,
    )
    if not np.allclose(coefficients, expected_coefficients, atol=5e-13):
        raise AssertionError(("ten theorem coefficients", coefficients))

    top = loads(str(audit.unresolved_orbits[0]["entries"]))
    if top[0] != [[3, 1, 1, 5], [1, 1, 0, 3]]:
        raise AssertionError(("top unproved physical orbit", top))
    if not np.isclose(audit.top_unresolved_gate, 0.0414623182965, atol=5e-13):
        raise AssertionError(("physical-orbit gate", audit.top_unresolved_gate))
    if not np.isclose(audit.top_provisional_gate, 0.0379251204234, atol=5e-13):
        raise AssertionError(("provisional-orbit gate", audit.top_provisional_gate))

    window = {int(row["N"]): row for row in audit.window}
    if window[512]["status"] != "unsupported_by_current_N_equals_q_squared_witness":
        raise AssertionError(("N=512 geometry", window[512]))
    if window[2048]["status"] != "unsupported_by_current_N_equals_q_squared_witness":
        raise AssertionError(("N=2048 geometry", window[2048]))
    if window[256]["q"] != 16 or window[4096]["q"] != 64:
        raise AssertionError(("square-order window", window))

    missing_sources = sorted(
        {
            str(row["source"])
            for row in audit.rows
            if not (ROOT / str(row["source"])).is_file()
        }
    )
    if missing_sources:
        raise AssertionError(("missing ledger provenance", missing_sources))

    artifact_names = (
        "finite_size_ledger_summary.json",
        "finite_size_window.csv",
        "n1024_balanced_ledger.csv",
        "n1024_unresolved_orbits.csv",
    )
    with TemporaryDirectory() as temporary_directory:
        generated = Path(temporary_directory)
        write_artifacts(audit, generated)
        mismatches = [
            name
            for name in artifact_names
            if (generated / name).read_bytes()
            != (ROOT / "artifacts" / name).read_bytes()
        ]
    if mismatches:
        raise AssertionError(("stale committed ledger artifacts", mismatches))

    print(
        "finite-size passive ledger passed: "
        f"balanced={audit.counts['balanced_entries']},"
        f"proved={audit.counts['proved_arbitrary_upper_entries']},"
        f"unresolved={audit.counts['unresolved_entries']},"
        f"physical_gate={audit.top_unresolved_gate:.12g},"
        f"provisional_gate={audit.top_provisional_gate:.12g},"
        f"total={audit.totals.ten_theorem_total:.12g},"
        f"margin={audit.totals.ten_theorem_margin:.12g},"
        f"artifacts={len(artifact_names)}"
    )


if __name__ == "__main__":
    main()
