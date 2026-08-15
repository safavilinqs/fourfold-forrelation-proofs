#!/usr/bin/env python3
"""Targeted q=8 masked screen for the three open quintic complement pairs.

This is diagnostic evidence, not an arbitrary-law upper bound.  It evaluates
complete deterministic translation-subspace orbit matrices with every
physical distinctness mask for one representative of each complement pair.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from masked_translation_subspace_screen import (  # noqa: E402
    artifact_text,
    diagnostic,
)


TEMPLATE_INDICES = (68, 70, 71)


def residual_screen():
    """Return the deterministic ten-trial q=8 residual screen."""

    return diagnostic(
        order=8,
        subspace_dimension=6,
        trials=10,
        selected_indices=TEMPLATE_INDICES,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = residual_screen()
    if arguments.output is not None:
        arguments.output.write_text(artifact_text(result), encoding="utf-8")
    riskiest = max(result.rows, key=lambda row: row.optimized_nuclear_lower)
    print(
        "q8 separated-quintic residual screen: "
        f"templates={result.templates_screened},"
        f"maximum_lower={result.maximum_optimized_lower:.12g},"
        f"maximum_orbit_upper={result.maximum_tangent_upper_on_orbit:.12g},"
        f"above_one={result.templates_above_one},"
        f"riskiest_split={riskiest.split},"
        "status=no_counterexample_diagnostic_only"
    )


if __name__ == "__main__":
    main()
