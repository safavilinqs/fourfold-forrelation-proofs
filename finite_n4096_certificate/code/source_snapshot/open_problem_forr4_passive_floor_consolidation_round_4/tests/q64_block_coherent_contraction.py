#!/usr/bin/env python3
"""Regression for the q64 block-coherent high-sector theorem."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_block_coherent_contraction import (  # noqa: E402
    artifact_text,
    block_coherent_coefficient,
    block_coherent_entries,
    diagnostic,
)
from q64_paper_target_gate import ORDER  # noqa: E402
import occupation_compatible_sector_optimization as occupation  # noqa: E402


def main() -> None:
    result = diagnostic()
    observed_discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.closed_entries,
        result.record_sector_bounds,
        result.maximum_entry_profile,
        result.maximum_entry_split,
    )
    expected_discrete = (
        64,
        4096,
        16_384,
        70,
        196,
        (3, 3, 3, 3),
        (3, 0, 3, 0),
    )
    if observed_discrete != expected_discrete:
        raise AssertionError(("block-coherent discrete result", observed_discrete))
    if result.minimum_coefficient_exact != "1/4096":
        raise AssertionError(("minimum exact coefficient", result))
    if result.maximum_coefficient_exact != "2609304163/39728800944":
        raise AssertionError(("maximum exact coefficient", result))

    observed = (
        result.minimum_coefficient,
        result.maximum_coefficient,
        result.previous_two_tier.total,
        result.block_coherent_inserted.total,
        result.block_coherent_inserted.beta,
        result.block_coherent_inserted.perron_upper,
        result.block_coherent_inserted.promise_loss,
        result.block_coherent_inserted.margin_to_one_third,
        result.margin_improvement,
    )
    expected = (
        0.000244140625,
        0.0656778986780387,
        0.3191811621612196,
        0.3094050070080503,
        0.7464869472876845,
        0.29368091490472226,
        0.015724092103328063,
        0.023928326325283023,
        0.00977615515316932,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("block-coherent numeric result", observed))

    entries = block_coherent_entries()
    coefficients = {entry: block_coherent_coefficient(*entry) for entry in entries}
    if len(coefficients) != 70 or not all(
        isinstance(value, Fraction) and value > 0 for value in coefficients.values()
    ):
        raise AssertionError("block-coherent exact coefficient map")
    cubic_target = occupation.endpoint_singleton_slice_energies(ORDER)[2] ** 0.5
    for (profile, _), coefficient in coefficients.items():
        target = cubic_target if 3 in profile else 0.5
        if float(coefficient) > target:
            raise AssertionError(
                ("coefficient misses q64 target", profile, coefficient)
            )

    committed = (ROOT / "artifacts" / "q64_block_coherent_contraction.json").read_text(
        encoding="utf-8"
    )
    if committed != artifact_text(result):
        raise AssertionError("stale q64 block-coherent artifact")

    print(
        "q64 block-coherent contraction passed: "
        f"entries={result.closed_entries},"
        f"record_sectors={result.record_sector_bounds},"
        f"coefficient_max={result.maximum_coefficient:.12g},"
        f"total={result.block_coherent_inserted.total:.12g},"
        f"margin={result.block_coherent_inserted.margin_to_one_third:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
