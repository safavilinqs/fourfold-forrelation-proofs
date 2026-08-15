#!/usr/bin/env python3
"""Regression for the inherited q64 chain-aware theorem insertion."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_chain_aware_insertion import (  # noqa: E402
    artifact_text,
    chain_aware_entries,
    diagnostic,
)
from q64_block_coherent_contraction import block_coherent_entries  # noqa: E402


def main() -> None:
    result = diagnostic()
    observed_discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.theorem_families,
        result.chain_aware_entries,
        result.overlap_with_block_coherent,
        result.total_proved_open_entries,
        result.remaining_open_entries,
        result.minimum_family,
        result.maximum_family,
    )
    expected_discrete = (
        64,
        4096,
        16_384,
        10,
        40,
        0,
        110,
        778,
        "adjacent_cubic_slice",
        "separated_endpoint_slice",
    )
    if observed_discrete != expected_discrete:
        raise AssertionError(("chain-aware discrete result", observed_discrete))

    observed = (
        result.minimum_coefficient,
        result.maximum_coefficient,
        result.block_coherent_only.total,
        result.chain_aware_inserted.total,
        result.chain_aware_inserted.beta,
        result.chain_aware_inserted.perron_upper,
        result.chain_aware_inserted.promise_loss,
        result.chain_aware_inserted.margin_to_one_third,
        result.margin_improvement,
    )
    expected = (
        0.00580989204377444,
        0.12397463639031407,
        0.3094050070080503,
        0.2960908671821436,
        0.7466670531429521,
        0.2810331927754359,
        0.015057674406707706,
        0.0372424661511897,
        0.013314139825906679,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("chain-aware numeric result", observed))

    entries = chain_aware_entries()
    if len(entries) != 40 or not all(len(entry[0]) == 4 for entry in entries):
        raise AssertionError("chain-aware coefficient map")
    overlap = set(entries).intersection(block_coherent_entries())
    if overlap:
        raise AssertionError(("block/chain theorem overlap", overlap))

    for family in result.families:
        source = (ROOT / family.source).resolve()
        if not source.is_file():
            raise AssertionError(("missing theorem provenance", family.name, source))

    committed = (ROOT / "artifacts" / "q64_chain_aware_insertion.json").read_text(
        encoding="utf-8"
    )
    if committed != artifact_text(result):
        raise AssertionError("stale q64 chain-aware insertion artifact")

    print(
        "q64 chain-aware insertion passed: "
        f"families={result.theorem_families},"
        f"new_entries={result.chain_aware_entries},"
        f"proved_entries={result.total_proved_open_entries},"
        f"total={result.chain_aware_inserted.total:.12g},"
        f"margin={result.chain_aware_inserted.margin_to_one_third:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
