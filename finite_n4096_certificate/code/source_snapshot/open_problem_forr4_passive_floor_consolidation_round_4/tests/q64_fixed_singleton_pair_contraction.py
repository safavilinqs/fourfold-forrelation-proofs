#!/usr/bin/env python3
"""Regression for the q64 fixed-singleton pair theorem."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_fixed_singleton_pair_contraction import (  # noqa: E402
    artifact_text,
    diagnostic,
    fixed_singleton_pair_coefficient,
    fixed_singleton_pair_entries,
    fixed_singleton_pair_row_energy,
)
from q64_reversed_middle_pair_contraction import (  # noqa: E402
    reversed_middle_pair_entries,
)


def main() -> None:
    result = diagnostic()
    observed_discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.closed_entries,
        result.cubic_pair_completions,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
    )
    expected_discrete = (64, 4096, 16_384, 4, 8_382_465, 232, 236, 652)
    if observed_discrete != expected_discrete:
        raise AssertionError(("fixed-singleton discrete result", observed_discrete))

    observed = (
        result.quintic_fixed_pair_energy,
        result.hadamard_squared_factor,
        result.universal_middle_maximum,
        result.row_energy_bound,
        result.coefficient,
        result.previous_routing.total,
        result.fixed_singleton_inserted.total,
        result.fixed_singleton_inserted.beta,
        result.fixed_singleton_inserted.perron_upper,
        result.fixed_singleton_inserted.promise_loss,
        result.fixed_singleton_inserted.margin_to_one_third,
        result.margin_improvement,
    )
    expected = (
        41.37997581845239,
        0.000244140625,
        0.0002640168970814132,
        0.005902900714455929,
        0.07683033720123796,
        0.32377678092061524,
        0.32303469500415644,
        0.7461938577336442,
        0.3061648799814117,
        0.01686981502274475,
        0.010298638329176879,
        0.0007420859164588012,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("fixed-singleton numeric result", observed))
    if not np.isclose(
        fixed_singleton_pair_coefficient() ** 2,
        fixed_singleton_pair_row_energy(),
        rtol=1e-13,
    ):
        raise AssertionError("fixed-singleton coefficient identity")

    entries = fixed_singleton_pair_entries()
    if len(entries) != 4 or set(entries).intersection(reversed_middle_pair_entries()):
        raise AssertionError("fixed-singleton entry partition")

    committed = (
        ROOT / "artifacts" / "q64_fixed_singleton_pair_contraction.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 fixed-singleton artifact")

    print(
        "q64 fixed-singleton pair contraction passed: "
        f"entries={result.closed_entries},"
        f"coefficient={result.coefficient:.12g},"
        f"total={result.fixed_singleton_inserted.total:.12g},"
        f"margin={result.fixed_singleton_inserted.margin_to_one_third:.12g},"
        f"margin_gain={result.margin_improvement:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
