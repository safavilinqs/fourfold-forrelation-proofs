#!/usr/bin/env python3
"""Regression for the q64 reversed middle-pair theorem."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_reversed_middle_pair_contraction import (  # noqa: E402
    artifact_text,
    diagnostic,
    reversed_middle_pair_coefficient,
    reversed_middle_pair_entries,
    reversed_middle_pair_row_energy,
)
from q64_shifted_middle_pair_contraction import (  # noqa: E402
    shifted_middle_pair_entries,
)


def main() -> None:
    result = diagnostic()
    observed_discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.closed_entries,
        result.cubic_completions,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
    )
    expected_discrete = (64, 4096, 16_384, 8, 4094, 224, 232, 656)
    if observed_discrete != expected_discrete:
        raise AssertionError(("reversed middle-pair discrete result", observed_discrete))

    observed = (
        result.quintic_fixed_pair_energy,
        result.record_one_middle_maximum,
        result.record_three_middle_maximum,
        result.row_energy_bound,
        result.coefficient,
        result.previous_routing.total,
        result.reversed_pair_inserted.total,
        result.reversed_pair_inserted.beta,
        result.reversed_pair_inserted.perron_upper,
        result.reversed_pair_inserted.promise_loss,
        result.reversed_pair_inserted.margin_to_one_third,
        result.margin_improvement,
    )
    expected = (
        41.37997581845239,
        0.0002640168970814132,
        1.180403414670998e-06,
        0.011808684408503778,
        0.10866777079016472,
        0.32504506334697836,
        0.32377678092061524,
        0.7461844487675892,
        0.30686888242313387,
        0.016907898497481356,
        0.009556552412718078,
        0.001268282426363121,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("reversed middle-pair numeric result", observed))
    if not np.isclose(
        reversed_middle_pair_coefficient() ** 2,
        reversed_middle_pair_row_energy(),
        rtol=1e-13,
    ):
        raise AssertionError("reversed middle-pair coefficient identity")

    entries = reversed_middle_pair_entries()
    if len(entries) != 8 or set(entries).intersection(shifted_middle_pair_entries()):
        raise AssertionError("reversed middle-pair entry partition")

    committed = (
        ROOT / "artifacts" / "q64_reversed_middle_pair_contraction.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 reversed middle-pair artifact")

    print(
        "q64 reversed middle-pair contraction passed: "
        f"entries={result.closed_entries},"
        f"coefficient={result.coefficient:.12g},"
        f"total={result.reversed_pair_inserted.total:.12g},"
        f"margin={result.reversed_pair_inserted.margin_to_one_third:.12g},"
        f"margin_gain={result.margin_improvement:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
