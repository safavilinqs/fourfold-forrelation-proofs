#!/usr/bin/env python3
"""Regression for the q64 decorated adjacent complete-row theorem."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_decorated_adjacent_row_insertion import (  # noqa: E402
    artifact_text,
    decorated_adjacent_coefficient,
    decorated_adjacent_entries,
    diagnostic,
    inserted_coefficients,
    local_complete_row_parameters,
    pre_decorated_quintic_entries,
    remaining_quintic_entries,
)
from q64_dual_endpoint_schur_insertion import (  # noqa: E402
    inserted_coefficients as previous_inserted_coefficients,
)


def orbit(entry):
    profile, split = entry
    complement = tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )
    reverse = tuple(reversed(profile))
    return frozenset(
        {
            (profile, split),
            (profile, complement),
            (reverse, tuple(reversed(split))),
            (reverse, tuple(reversed(complement))),
        }
    )


def main() -> None:
    result = diagnostic()
    discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.closed_entries,
        result.closed_orbits,
        result.degree_twelve_entries,
        result.extreme_entries,
        result.balanced_entries,
        result.fixed_one_cubic_fixed_four_quintic_entries,
        result.fixed_two_cubic_fixed_three_quintic_entries,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
        result.remaining_quintic_entries,
        result.remaining_extreme_entries,
        result.remaining_balanced_entries,
    )
    expected_discrete = (
        64,
        4096,
        16_384,
        16,
        4,
        16,
        8,
        8,
        8,
        8,
        276,
        292,
        596,
        136,
        96,
        40,
    )
    if discrete != expected_discrete:
        raise AssertionError(("decorated-row discrete result", discrete))

    observed = (
        result.adjacent_row_energy_bound,
        result.coefficient,
        result.previous_routing.total,
        result.decorated_row_inserted.total,
        result.decorated_row_inserted.beta,
        result.decorated_row_inserted.perron_upper,
        result.decorated_row_inserted.promise_loss,
        result.decorated_row_inserted.margin_to_one_third,
        result.routing_margin_improvement,
        result.remaining_quintic_local_proxy.total,
        result.remaining_quintic_local_proxy.margin_to_one_third,
        result.proxy_reserve_after_declared_allowance,
    )
    expected = (
        0.0004031890208244533,
        0.020079567246941685,
        0.32988360597454097,
        0.32910720773187735,
        0.7460989053986221,
        0.3118492405250244,
        0.017257967206852944,
        0.004226125601455966,
        0.0007763982426636229,
        0.33007665086528165,
        0.0032566824680516637,
        0.0022566824680516637,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("decorated-row numeric result", observed))
    if not np.isclose(
        decorated_adjacent_coefficient() ** 2,
        result.adjacent_row_energy_bound,
        rtol=1e-13,
    ):
        raise AssertionError("decorated-row coefficient identity")

    entries = decorated_adjacent_entries()
    expected_generators = {
        ((1, 3, 5, 3), (0, 1, 2, 3)),
        ((1, 3, 5, 3), (0, 2, 1, 3)),
        ((3, 1, 3, 5), (0, 1, 1, 4)),
        ((3, 1, 3, 5), (0, 1, 2, 3)),
    }
    observed_orbits = {orbit(entry) for entry in entries}
    if observed_orbits != {orbit(entry) for entry in expected_generators}:
        raise AssertionError(("decorated-row orbit inventory", observed_orbits))
    if any(local_complete_row_parameters(entry) is None for entry in entries):
        raise AssertionError("decorated-row topology predicate")
    for entry in entries:
        profile, split = entry
        start = local_complete_row_parameters(entry)[0]  # type: ignore[index]
        outer = 3 if start == 0 else 0
        if profile[outer] != 3 or split[outer] not in (0, 3):
            raise AssertionError(("decorated row has split outer block", entry))
    parameters = {
        local_complete_row_parameters(entry)[1:]  # type: ignore[index]
        for entry in entries
    }
    if parameters != {(1, 4), (2, 3)}:
        raise AssertionError(("decorated-row fixed counts", parameters))

    pre_decorated = set(pre_decorated_quintic_entries())
    remaining = set(remaining_quintic_entries())
    if not set(entries).issubset(pre_decorated):
        raise AssertionError("decorated-row entries outside input family")
    if set(entries).intersection(remaining):
        raise AssertionError("decorated-row partition overlap")
    if len(remaining) + len(entries) != len(pre_decorated):
        raise AssertionError("decorated-row partition size")
    for entry in entries:
        if np.isclose(
            inserted_coefficients()[entry],
            previous_inserted_coefficients()[entry],
        ):
            raise AssertionError(("decorated-row theorem not inserted", entry))

    committed = (
        ROOT / "artifacts" / "q64_decorated_adjacent_row_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 decorated-row artifact")

    print(
        "q64 decorated adjacent-row insertion passed: "
        f"entries={result.closed_entries},"
        f"coefficient={result.coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.decorated_row_inserted.total:.12g},"
        f"margin={result.decorated_row_inserted.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
