#!/usr/bin/env python3
"""Regression for the q=2 arbitrary-diagonal alternating certificate."""

from __future__ import annotations

from pathlib import Path
import sys


SEARCHES = Path(__file__).resolve().parents[1] / "searches"
sys.path.insert(0, str(SEARCHES))

from weighted_alternating_q2_certificate import certificate


def main() -> None:
    result = certificate()
    if result.supporting_upper >= 0.471845:
        raise AssertionError(("weighted q2 upper", result.supporting_upper))
    if result.high_columns != 48 or result.low_columns != 96:
        raise AssertionError(
            ("weighted q2 orbit sizes", result.high_columns, result.low_columns)
        )
    print(
        "weighted alternating q=2 regression passed: "
        f"uniform={result.uniform_value:.12g},"
        f"optimum={result.orbit_value},"
        f"upper={result.supporting_upper:.12g}"
    )


if __name__ == "__main__":
    main()
