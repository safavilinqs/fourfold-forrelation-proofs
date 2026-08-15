#!/usr/bin/env python3
"""Regression for the q64 whole-cubic decorated completion rows."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import numpy as np
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_whole_cubic_decorated_row_insertion import (  # noqa: E402
    CUBIC_LEADING_CUBIC_QUINTIC,
    CUBIC_LEADING_QUINTIC_CUBIC,
    CUBIC_QUINTIC_TRAILING_CUBIC,
    QUINTIC_CUBIC_TRAILING_CUBIC,
    artifact_text,
    coefficient_map,
    cubic_endpoint_energy_parts,
    cubic_leading_coefficient,
    cubic_trailing_coefficient,
    diagnostic,
    inserted_coefficients,
    pre_whole_cubic_quintic_entries,
    quintic_endpoint_energy_parts,
    quintic_leading_coefficient,
    quintic_trailing_coefficient,
    remaining_quintic_entries,
    sector_parameters,
    whole_cubic_decorated_entries,
)
from q64_degree_ten_completion_row_insertion import (  # noqa: E402
    inserted_coefficients as previous_inserted_coefficients,
    orbit,
)
from adjacent_cubic_quintic_orbit_witness import (  # noqa: E402
    exact_link_moments,
    parity_record_size,
)
from opposite_endpoint_orbit_scan import (  # noqa: E402
    cubic_weight,
    quintic_weight,
)


def transpose_support(
    order: int,
    support: tuple[int, ...],
) -> tuple[int, ...]:
    return tuple(
        sorted((cell % order) * order + cell // order for cell in support)
    )


def q4_scalar_row_maxima() -> tuple[tuple[float, float, float], ...]:
    q = 4
    dimension = q * q
    moments = exact_link_moments(q)
    cubics = moments.supports_three
    quintics = moments.supports_five
    pairs = tuple(combinations(range(dimension), 2))
    fours = tuple(combinations(range(dimension), 4))
    pair_index = {support: index for index, support in enumerate(pairs)}
    four_index = {support: index for index, support in enumerate(fours)}

    cubic_pair = sparse.lil_matrix((len(pairs), len(cubics)))
    for position, support in enumerate(cubics):
        for fixed in combinations(support, 2):
            cubic_pair[pair_index[fixed], position] = 1
    quintic_four = sparse.lil_matrix((len(quintics), len(fours)))
    for position, support in enumerate(quintics):
        for fixed in combinations(support, 4):
            quintic_four[position, four_index[fixed]] = 1

    cubic_endpoint = np.asarray(
        [
            cubic_weight(transpose_support(q, support), q)
            for support in cubics
        ]
    )
    quintic_endpoint = np.asarray(
        [quintic_weight(support, q) for support in quintics]
    )
    records = np.asarray(
        [parity_record_size(q, support, axis=1) for support in cubics]
    )
    moment_squared = np.square(moments.moment_35)

    cubic_tables = []
    quintic_tables = []
    for record in (1, 3):
        sector = (records == record)[:, None]
        cubic_weighted = (
            moment_squared * np.square(cubic_endpoint)[:, None] * sector
        )
        quintic_weighted = (
            moment_squared * np.square(quintic_endpoint)[None, :] * sector
        )
        cubic_tables.append(
            np.asarray(
                cubic_pair.tocsr()
                @ cubic_weighted
                @ quintic_four.tocsr()
            )
        )
        quintic_tables.append(
            np.asarray(
                cubic_pair.tocsr()
                @ quintic_weighted
                @ quintic_four.tocsr()
            )
        )

    def maxima(tables: list[np.ndarray]) -> tuple[float, float, float]:
        return (
            float(tables[0].max()),
            float(tables[1].max()),
            float((tables[0] + tables[1]).max()),
        )

    return maxima(cubic_tables), maxima(quintic_tables)


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
        304,
        320,
        568,
        108,
        72,
        36,
    )
    if discrete != expected_discrete:
        raise AssertionError(("whole-cubic discrete result", discrete))

    observed = (
        result.record_one_middle_maximum,
        result.record_three_middle_maximum,
        result.cubic_endpoint_record_one_energy,
        result.cubic_endpoint_record_three_energy,
        result.cubic_leading_residual_factor,
        result.cubic_leading_coefficient,
        result.cubic_trailing_residual_factor,
        result.cubic_trailing_coefficient,
        result.quintic_endpoint_record_one_energy,
        result.quintic_endpoint_record_three_energy,
        result.quintic_leading_residual_factor,
        result.quintic_leading_coefficient,
        result.quintic_trailing_residual_factor,
        result.quintic_trailing_coefficient,
        result.previous_routing.total,
        result.decorated_rows_inserted.total,
        result.decorated_rows_inserted.beta,
        result.decorated_rows_inserted.perron_upper,
        result.decorated_rows_inserted.promise_loss,
        result.decorated_rows_inserted.margin_to_one_third,
        result.routing_change,
        result.remaining_quintic_local_proxy.total,
        result.remaining_quintic_local_proxy.margin_to_one_third,
        result.proxy_reserve_after_declared_allowance,
    )
    expected = (
        0.0002640168970814132,
        2.4001536098310292e-05,
        9.055001288789966e-06,
        0.0001485469219689924,
        0.015625,
        0.00019615563220401556,
        1.0,
        0.012553960461056996,
        0.03593930011520738,
        0.00965076050581523,
        0.015625,
        0.0033362232979529936,
        1.0,
        0.2135182910689916,
        0.32528397960750105,
        0.32528397960750105,
        0.7461469557307886,
        0.3082235000172111,
        0.017060479590289935,
        0.008049353725832264,
        0.0,
        0.32598862373932996,
        0.007344709594003351,
        0.006344709594003351,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("whole-cubic numeric result", observed))

    coefficient_identities = (
        (cubic_leading_coefficient(), cubic_endpoint_energy_parts(), 1 / 64),
        (cubic_trailing_coefficient(), cubic_endpoint_energy_parts(), 1),
        (
            quintic_leading_coefficient(),
            quintic_endpoint_energy_parts(),
            1 / 64,
        ),
        (quintic_trailing_coefficient(), quintic_endpoint_energy_parts(), 1),
    )
    for coefficient, parts, residual in coefficient_identities:
        if not np.isclose(
            coefficient**2,
            sum(parts) * residual**2,
            rtol=1e-13,
        ):
            raise AssertionError(("whole-cubic coefficient identity", coefficient))
    if sector_parameters() != (
        result.record_one_middle_maximum,
        result.record_three_middle_maximum,
    ):
        raise AssertionError("whole-cubic sector parameters")

    generators = (
        CUBIC_LEADING_CUBIC_QUINTIC,
        CUBIC_QUINTIC_TRAILING_CUBIC,
        CUBIC_LEADING_QUINTIC_CUBIC,
        QUINTIC_CUBIC_TRAILING_CUBIC,
    )
    expected_entries = set().union(*(set(orbit(entry)) for entry in generators))
    entries = set(whole_cubic_decorated_entries())
    if entries != expected_entries or len(entries) != 16:
        raise AssertionError(("whole-cubic orbit inventory", entries))
    pre_whole = set(pre_whole_cubic_quintic_entries())
    remaining = set(remaining_quintic_entries())
    if not entries.issubset(pre_whole):
        raise AssertionError("whole-cubic entries outside live family")
    if entries.intersection(remaining):
        raise AssertionError("whole-cubic partition overlap")
    if len(entries) + len(remaining) != len(pre_whole):
        raise AssertionError("whole-cubic partition size")
    if set(coefficient_map()) != entries:
        raise AssertionError("whole-cubic coefficient map")
    for entry in entries:
        if np.isclose(
            inserted_coefficients()[entry],
            previous_inserted_coefficients()[entry],
        ) and not np.isclose(
            coefficient_map()[entry],
            previous_inserted_coefficients()[entry],
        ):
            raise AssertionError(("whole-cubic theorem not inserted", entry))

    q4 = q4_scalar_row_maxima()
    expected_q4 = (
        (0.20370370370370364, 1.5555555555555562, 1.5555555555555562),
        (0.4259259259259259, 2.0000000000000004, 2.0000000000000004),
    )
    if not np.allclose(q4, expected_q4, rtol=1e-13, atol=1e-14):
        raise AssertionError(("whole-cubic q4 exact rows", q4))

    committed = (
        ROOT / "artifacts" / "q64_whole_cubic_decorated_row_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 whole-cubic decorated-row artifact")

    print(
        "q64 whole-cubic decorated-row insertion passed: "
        f"entries={result.closed_entries},"
        f"coefficients={result.cubic_leading_coefficient:.12g}/"
        f"{result.cubic_trailing_coefficient:.12g}/"
        f"{result.quintic_leading_coefficient:.12g}/"
        f"{result.quintic_trailing_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.decorated_rows_inserted.total:.12g},"
        f"margin={result.decorated_rows_inserted.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
