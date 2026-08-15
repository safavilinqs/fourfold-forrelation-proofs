#!/usr/bin/env python3
"""Regression for the current degree-six joint-occupation obstruction."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from degree_six_joint_occupation_optimization import certificate


def main() -> None:
    result = certificate()
    available_margin = 0.160358131958
    if not result.objective <= result.supporting_upper:
        raise AssertionError(("supporting upper", result))
    if not 0.349 <= result.supporting_upper <= 0.351:
        raise AssertionError(("joint upper", result))
    if not 0.245 <= result.degree_six_at_candidate <= 0.247:
        raise AssertionError(("degree-six attribution", result))
    if not 0.103 <= result.double_endpoint_at_candidate <= 0.105:
        raise AssertionError(("double endpoint attribution", result))
    if not 0.188 <= result.supporting_upper - available_margin <= 0.191:
        raise AssertionError(("barrier size", result))
    print(
        "degree-six joint occupation barrier confirmed: "
        f"supporting_upper={result.supporting_upper:.12g},"
        f"degree_six={result.degree_six_at_candidate:.12g},"
        f"double_endpoint={result.double_endpoint_at_candidate:.12g},"
        f"overshoot={result.supporting_upper-available_margin:.12g}"
    )


if __name__ == "__main__":
    main()
