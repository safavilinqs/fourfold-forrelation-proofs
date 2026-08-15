#!/usr/bin/env python3
"""Regression for the q64 shifted middle-pair theorem."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_chain_aware_insertion import chain_aware_entries  # noqa: E402
from q64_shifted_middle_pair_contraction import (  # noqa: E402
    artifact_text,
    diagnostic,
    shifted_middle_pair_coefficient,
    shifted_middle_pair_entries,
    shifted_middle_pair_row_energy,
)


def main() -> None:
    result = diagnostic()
    observed_discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.closed_entries,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
    )
    expected_discrete = (64, 4096, 16_384, 4, 220, 224, 664)
    if observed_discrete != expected_discrete:
        raise AssertionError(("shifted middle-pair discrete result", observed_discrete))

    observed = (
        result.cubic_fixed_singleton_energy,
        result.quintic_fixed_triple_energy,
        result.record_one_middle_maximum,
        result.record_three_middle_maximum,
        result.row_energy_bound,
        result.coefficient,
        result.previous_routing.total,
        result.shifted_pair_inserted.total,
        result.shifted_pair_inserted.beta,
        result.shifted_pair_inserted.perron_upper,
        result.shifted_pair_inserted.promise_loss,
        result.shifted_pair_inserted.margin_to_one_third,
        result.margin_improvement,
    )
    expected = (
        0.500244140625,
        1.4538457961309523,
        0.0002640168970814132,
        1.180403414670998e-06,
        0.00020764608565648136,
        0.014409930105884669,
        0.33193582943438404,
        0.32504506334697836,
        0.7461670784772364,
        0.3080666393506051,
        0.016978423996373303,
        0.008288269986354957,
        0.006890766087405686,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("shifted middle-pair numeric result", observed))
    if not np.isclose(
        shifted_middle_pair_coefficient() ** 2,
        shifted_middle_pair_row_energy(),
        rtol=1e-13,
    ):
        raise AssertionError("shifted middle-pair coefficient identity")

    entries = shifted_middle_pair_entries()
    if len(entries) != 4 or set(entries).intersection(chain_aware_entries()):
        raise AssertionError("shifted middle-pair entry partition")

    committed = (
        ROOT / "artifacts" / "q64_shifted_middle_pair_contraction.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 shifted middle-pair artifact")

    print(
        "q64 shifted middle-pair contraction passed: "
        f"entries={result.closed_entries},"
        f"coefficient={result.coefficient:.12g},"
        f"total={result.shifted_pair_inserted.total:.12g},"
        f"margin={result.shifted_pair_inserted.margin_to_one_third:.12g},"
        f"margin_gain={result.margin_improvement:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
