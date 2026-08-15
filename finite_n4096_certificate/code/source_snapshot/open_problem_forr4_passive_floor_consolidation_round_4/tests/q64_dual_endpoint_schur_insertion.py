#!/usr/bin/env python3
"""Regression for the q64 dual endpoint-slice Schur theorem."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_dual_endpoint_schur_insertion import (  # noqa: E402
    artifact_text,
    diagnostic,
    dual_endpoint_entries,
    has_favorable_cubic_singleton,
    has_favorable_quintic_singleton,
    inserted_coefficients,
    local_slice_coefficients,
    pre_dual_quintic_entries,
    remaining_quintic_entries,
    split_cubic_index,
)
from q64_fixed_pair_adjacent_row_contraction import (  # noqa: E402
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
        result.degree_ten_entries,
        result.degree_twelve_entries,
        result.extreme_entries,
        result.balanced_entries,
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
        12,
        3,
        4,
        8,
        0,
        12,
        264,
        276,
        612,
        152,
        104,
        48,
    )
    if discrete != expected_discrete:
        raise AssertionError(("dual-endpoint discrete result", discrete))

    observed = (
        result.cubic_fixed_pair_energy,
        result.cubic_schur_factor,
        result.balanced_quintic_fixed_triple_energy,
        result.balanced_quintic_schur_factor,
        result.balanced_coefficient,
        result.previous_routing.total,
        result.dual_endpoint_inserted.total,
        result.dual_endpoint_inserted.beta,
        result.dual_endpoint_inserted.perron_upper,
        result.dual_endpoint_inserted.promise_loss,
        result.dual_endpoint_inserted.margin_to_one_third,
        result.routing_margin_spent,
        result.remaining_quintic_local_proxy.total,
        result.remaining_quintic_local_proxy.margin_to_one_third,
        result.proxy_reserve_after_declared_allowance,
    )
    expected = (
        0.015384734623015874,
        0.12403521525363623,
        1.4538457961309523,
        1.2057552803661913,
        0.14955611574342903,
        0.3292482541342017,
        0.32988360597454097,
        0.7460871349070535,
        0.3125769266646106,
        0.01730667930993033,
        0.0034497273587923427,
        0.000635351840339249,
        0.33106696650853623,
        0.0022663668247970836,
        0.0012663668247970836,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("dual-endpoint numeric result", observed))
    if not np.isclose(
        result.cubic_schur_factor
        * result.balanced_quintic_schur_factor,
        result.balanced_coefficient,
        rtol=1e-13,
    ):
        raise AssertionError("dual-endpoint coefficient product")
    if local_slice_coefficients()[1] != result.balanced_coefficient:
        raise AssertionError("dual-endpoint local coefficient")

    entries = dual_endpoint_entries()
    if len({orbit(entry) for entry in entries}) != 3:
        raise AssertionError("dual-endpoint orbit count")
    previous = set(previous_inserted_coefficients())
    if not all(
        has_favorable_cubic_singleton(entry)
        and has_favorable_quintic_singleton(entry)
        and split_cubic_index(entry) in range(4)
        for entry in entries
    ):
        raise AssertionError("dual-endpoint topology predicate")
    pre_dual = set(pre_dual_quintic_entries())
    remaining = set(remaining_quintic_entries())
    if not set(entries).issubset(pre_dual) or set(entries).intersection(remaining):
        raise AssertionError("dual-endpoint partition overlap")
    if len(remaining) + len(entries) != len(pre_dual):
        raise AssertionError("dual-endpoint partition size")
    for entry in entries:
        if np.isclose(
            inserted_coefficients()[entry],
            previous_inserted_coefficients()[entry],
        ):
            raise AssertionError(("dual-endpoint theorem not inserted", entry))
    if len(previous) != 6016:
        raise AssertionError(("unexpected q64 coefficient map", len(previous)))

    committed = (
        ROOT / "artifacts" / "q64_dual_endpoint_schur_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 dual-endpoint artifact")

    print(
        "q64 dual-endpoint Schur insertion passed: "
        f"entries={result.closed_entries},"
        f"coefficient={result.balanced_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.dual_endpoint_inserted.total:.12g},"
        f"margin={result.dual_endpoint_inserted.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
