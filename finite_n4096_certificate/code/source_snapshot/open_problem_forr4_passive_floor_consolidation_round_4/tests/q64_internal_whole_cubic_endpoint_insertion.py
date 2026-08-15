#!/usr/bin/env python3
"""Regression for the q64 internal whole-cubic endpoint theorem."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_internal_whole_cubic_endpoint_insertion import (  # noqa: E402
    artifact_text,
    coefficients,
    diagnostic,
    distinctness_factors,
    inserted_coefficients,
    internal_whole_cubic_entries,
    pre_internal_whole_quintic_entries,
    remaining_quintic_entries,
    same_side_adjacent_whole_cubic,
)
from q64_last_degree_ten_chain_insertion import (  # noqa: E402
    inserted_coefficients as previous_inserted_coefficients,
)
from q64_degree_ten_completion_row_insertion import orbit  # noqa: E402
from opposite_endpoint_orbit_scan import endpoint_moment  # noqa: E402


EXPECTED_GENERATORS = (
    ((1, 3, 3, 5), (0, 0, 2, 4)),
    ((1, 3, 5, 3), (0, 0, 4, 2)),
    ((3, 1, 3, 5), (1, 1, 3, 1)),
    ((3, 3, 1, 5), (1, 3, 1, 1)),
)


def endpoint_magnitude_check(order: int) -> int:
    dimension = order * order
    checked = 0
    supports = combinations(range(dimension), 3)
    for support in supports:
        for singleton in range(0, dimension, max(1, dimension // 7)):
            moment = endpoint_moment(
                support,
                singleton,
                order,
                3,
                False,
            )
            if abs(moment) > 1 / order + 2e-14:
                raise AssertionError(
                    ("whole-cubic endpoint magnitude", order, support, moment)
                )
            checked += 1
        if checked >= 10_000:
            break
    return checked


def main() -> None:
    result = diagnostic()
    discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.closed_entries,
        result.closed_orbits,
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
        16,
        4,
        16,
        0,
        324,
        340,
        548,
        88,
        56,
        32,
    )
    if discrete != expected_discrete:
        raise AssertionError(("internal whole-cubic discrete result", discrete))

    observed = (
        result.cubic_distinctness_factor,
        result.extreme_quintic_distinctness_factor,
        result.balanced_quintic_distinctness_factor,
        result.endpoint_factor,
        result.extreme_coefficient,
        result.balanced_coefficient,
        result.previous_routing.total,
        result.endpoint_inserted.total,
        result.endpoint_inserted.beta,
        result.endpoint_inserted.perron_upper,
        result.endpoint_inserted.promise_loss,
        result.endpoint_inserted.margin_to_one_third,
        result.routing_change,
        result.remaining_quintic_local_proxy.total,
        result.remaining_quintic_local_proxy.margin_to_one_third,
        result.proxy_reserve_after_declared_allowance,
    )
    expected = (
        2.413207370108021,
        2.9978026001826406,
        7.455480652750664,
        0.015625,
        0.11303623951390214,
        0.2811190759205572,
        0.3247663269533303,
        0.32476632695333046,
        0.7461529246324944,
        0.3077302268922124,
        0.017036100061118036,
        0.00856700638000285,
        1.6653345369377348e-16,
        0.3250376129305583,
        0.008295720402775042,
        0.007295720402775042,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("internal whole-cubic numeric result", observed))
    factors = distinctness_factors()
    coefficient_pair = coefficients()
    if not np.allclose(
        coefficient_pair,
        (factors[0] * factors[1] / 64, factors[0] * factors[2] / 64),
        rtol=1e-13,
    ):
        raise AssertionError("internal whole-cubic coefficient identity")

    expected_entries = set().union(
        *(set(orbit(entry)) for entry in EXPECTED_GENERATORS)
    )
    entries = set(internal_whole_cubic_entries())
    if entries != expected_entries or len(entries) != 16:
        raise AssertionError(("internal whole-cubic orbit inventory", entries))
    if not all(same_side_adjacent_whole_cubic(entry) for entry in entries):
        raise AssertionError("internal whole-cubic topology predicate")
    pre_internal = set(pre_internal_whole_quintic_entries())
    remaining = set(remaining_quintic_entries())
    if not entries.issubset(pre_internal) or entries.intersection(remaining):
        raise AssertionError("internal whole-cubic partition overlap")
    if len(entries) + len(remaining) != len(pre_internal):
        raise AssertionError("internal whole-cubic partition size")
    for entry in entries:
        if np.isclose(
            inserted_coefficients()[entry],
            previous_inserted_coefficients()[entry],
        ):
            raise AssertionError(("internal whole-cubic theorem not inserted", entry))

    endpoint_checks = sum(endpoint_magnitude_check(order) for order in (4, 8))

    committed = (
        ROOT
        / "artifacts"
        / "q64_internal_whole_cubic_endpoint_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 internal whole-cubic endpoint artifact")

    print(
        "q64 internal whole-cubic endpoint insertion passed: "
        f"entries={result.closed_entries},"
        f"coefficients={result.extreme_coefficient:.12g}/"
        f"{result.balanced_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.endpoint_inserted.total:.12g},"
        f"margin={result.endpoint_inserted.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        f"endpoint_checks={endpoint_checks},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
