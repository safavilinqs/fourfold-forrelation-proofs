#!/usr/bin/env python3
"""Regression for the full 64-cut double-endpoint occupation ledger."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from double_endpoint_occupation_optimization import certificate, coefficient


def main() -> None:
    expected_whole_block = {
        (0, 0, 0, 0): 32**-3,
        (3, 0, 0, 0): 32**-2,
        (0, 1, 0, 0): 32**-2,
        (3, 1, 0, 0): 32**-2,
        (3, 0, 1, 0): 32**-1,
        (0, 1, 1, 0): 32**-1,
    }
    for split, expected in expected_whole_block.items():
        if coefficient(split) != expected:
            raise AssertionError(("whole-block cut table", split))
    result = certificate()
    available_margin = 0.160358131958
    if not result.objective <= result.supporting_upper:
        raise AssertionError(("supporting upper", result))
    if not 0.115 < result.attenuated_upper < 0.117:
        raise AssertionError(("attenuated full ledger", result))
    if result.attenuated_upper >= available_margin:
        raise AssertionError(("N=1024 double-endpoint gate", result))
    print(
        "double-endpoint occupation optimization passed: "
        f"objective={result.objective:.12g},"
        f"supporting_upper={result.supporting_upper:.12g},"
        f"attenuated_upper={result.attenuated_upper:.12g},"
        f"margin_slack={available_margin-result.attenuated_upper:.12g}"
    )


if __name__ == "__main__":
    main()
