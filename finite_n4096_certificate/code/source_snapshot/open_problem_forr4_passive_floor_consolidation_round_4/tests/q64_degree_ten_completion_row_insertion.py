#!/usr/bin/env python3
"""Regression for the q64 degree-ten completion-row theorems."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import numpy as np
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_degree_ten_completion_row_insertion import (  # noqa: E402
    DOUBLE_ENDPOINT_ONE_FOUR,
    LEFT_DOUBLE_SINGLETON,
    REVERSED_DOUBLE_SINGLETON,
    artifact_text,
    coefficient_map,
    degree_ten_completion_entries,
    diagnostic,
    double_endpoint_one_four_coefficient,
    double_endpoint_one_four_energy_parts,
    inserted_coefficients,
    left_double_singleton_coefficient,
    left_double_singleton_energy_parts,
    orbit,
    pre_degree_ten_quintic_entries,
    remaining_quintic_entries,
    reversed_double_singleton_coefficient,
    reversed_double_singleton_energy_parts,
    sector_parameters,
)
from q64_decorated_adjacent_row_insertion import (  # noqa: E402
    inserted_coefficients as previous_inserted_coefficients,
)
from adjacent_cubic_quintic_orbit_witness import (  # noqa: E402
    exact_link_moments,
    parity_record_size,
)
from opposite_endpoint_orbit_scan import (  # noqa: E402
    cubic_weight,
    quintic_weight,
    sylvester,
)


def transpose_support(
    order: int,
    support: tuple[int, ...],
) -> tuple[int, ...]:
    return tuple(
        sorted((cell % order) * order + cell // order for cell in support)
    )


def q4_scalar_row_maxima() -> tuple[tuple[float, float, float], ...]:
    """Compute exact record-one, record-three, and total row maxima."""

    q = 4
    dimension = q * q
    moments = exact_link_moments(q)
    cubics = moments.supports_three
    quintics = moments.supports_five
    pairs = tuple(combinations(range(dimension), 2))
    triples = cubics
    fours = tuple(combinations(range(dimension), 4))
    pair_index = {support: index for index, support in enumerate(pairs)}
    triple_index = {support: index for index, support in enumerate(triples)}
    four_index = {support: index for index, support in enumerate(fours)}

    cubic_singleton = sparse.lil_matrix((dimension, len(cubics)))
    cubic_pair = sparse.lil_matrix((len(pairs), len(cubics)))
    for position, support in enumerate(cubics):
        for cell in support:
            cubic_singleton[cell, position] = 1
        for fixed in combinations(support, 2):
            cubic_pair[pair_index[fixed], position] = 1

    quintic_triple = sparse.lil_matrix((len(quintics), len(triples)))
    quintic_four = sparse.lil_matrix((len(quintics), len(fours)))
    for position, support in enumerate(quintics):
        for fixed in combinations(support, 3):
            quintic_triple[position, triple_index[fixed]] = 1
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

    left_tables = []
    reversed_tables = []
    endpoint_tables = []
    for record in (1, 3):
        sector = (records == record)[:, None]
        left_weighted = (
            moment_squared * np.square(cubic_endpoint)[:, None] * sector
        )
        reversed_weighted = (
            moment_squared * np.square(quintic_endpoint)[None, :] * sector
        )
        endpoint_weighted = (
            left_weighted * np.square(quintic_endpoint)[None, :]
        )
        left_tables.append(
            np.asarray(
                cubic_pair.tocsr()
                @ left_weighted
                @ quintic_triple.tocsr()
            )
        )
        reversed_tables.append(
            np.asarray(
                cubic_singleton.tocsr()
                @ reversed_weighted
                @ quintic_four.tocsr()
            )
        )
        endpoint_tables.append(
            np.asarray(
                cubic_singleton.tocsr()
                @ endpoint_weighted
                @ quintic_four.tocsr()
            )
        )

    def maxima(tables: list[np.ndarray]) -> tuple[float, float, float]:
        return (
            float(tables[0].max()),
            float(tables[1].max()),
            float((tables[0] + tables[1]).max()),
        )

    return maxima(left_tables), maxima(reversed_tables), maxima(endpoint_tables)


def double_singleton_residual_stress() -> float:
    """Stress the claimed ``1/q`` Walsh residual under arbitrary laws."""

    rng = np.random.default_rng(2026071611)
    q = 4
    dimension = q * q
    hadamard = sylvester(dimension)
    worst = 0.0
    for _ in range(32):
        rows = rng.integers(dimension, size=28)
        columns = rng.integers(dimension, size=(44, 3))
        matrix = np.empty((len(rows), len(columns)))
        for row_index, fixed_xor in enumerate(rows):
            for column_index, (first, middle, twist) in enumerate(columns):
                matrix[row_index, column_index] = (
                    hadamard[first, middle]
                    * hadamard[fixed_xor ^ twist, middle]
                )
        row_law = rng.dirichlet(np.ones(len(rows)))
        column_law = rng.dirichlet(np.ones(len(columns)))
        weighted = (
            np.sqrt(row_law)[:, None]
            * matrix
            * np.sqrt(column_law)[None, :]
        )
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > 1 / q + 5e-12:
            raise AssertionError(("double-singleton Walsh residual", nuclear))
        worst = max(worst, nuclear)
    return worst


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
        12,
        3,
        8,
        4,
        292,
        304,
        584,
        124,
        88,
        36,
    )
    if discrete != expected_discrete:
        raise AssertionError(("degree-ten discrete result", discrete))

    observed = (
        result.record_one_middle_maximum,
        result.record_three_middle_maximum,
        result.left_double_singleton_record_one_energy,
        result.left_double_singleton_record_three_energy,
        result.left_double_singleton_residual_factor,
        result.left_double_singleton_coefficient,
        result.reversed_double_singleton_record_one_energy,
        result.reversed_double_singleton_record_three_energy,
        result.reversed_double_singleton_residual_factor,
        result.reversed_double_singleton_coefficient,
        result.double_endpoint_record_one_energy,
        result.double_endpoint_record_three_energy,
        result.double_endpoint_coefficient,
        result.previous_routing.total,
        result.degree_ten_inserted.total,
        result.degree_ten_inserted.beta,
        result.degree_ten_inserted.perron_upper,
        result.degree_ten_inserted.promise_loss,
        result.degree_ten_inserted.margin_to_one_third,
        result.routing_margin_improvement,
        result.remaining_quintic_local_proxy.total,
        result.remaining_quintic_local_proxy.margin_to_one_third,
        result.proxy_reserve_after_declared_allowance,
    )
    expected = (
        0.0002640168970814132,
        2.4001536098310292e-05,
        3.484803526291896e-05,
        0.304001275809543,
        0.015625,
        0.008615542310153686,
        3.396263860887097,
        19.759932135656683,
        0.015625,
        0.07518888324226325,
        0.0008556976217906518,
        0.004830094840355497,
        0.07540419392942377,
        0.32910720773187735,
        0.32528397960750105,
        0.7461469557309,
        0.30822350001766624,
        0.01706047958983481,
        0.008049353725832264,
        0.003823228124376299,
        0.32598862373933,
        0.007344709594003296,
        0.006344709594003296,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("degree-ten numeric result", observed))

    coefficient_identities = (
        (
            left_double_singleton_coefficient(),
            left_double_singleton_energy_parts(),
            1 / 64,
        ),
        (
            reversed_double_singleton_coefficient(),
            reversed_double_singleton_energy_parts(),
            1 / 64,
        ),
        (
            double_endpoint_one_four_coefficient(),
            double_endpoint_one_four_energy_parts(),
            1,
        ),
    )
    for coefficient, parts, residual in coefficient_identities:
        if not np.isclose(
            coefficient**2,
            sum(parts) * residual**2,
            rtol=1e-13,
        ):
            raise AssertionError(("degree-ten coefficient identity", coefficient))
    if sector_parameters() != (
        result.record_one_middle_maximum,
        result.record_three_middle_maximum,
    ):
        raise AssertionError("degree-ten sector parameters")

    expected_entries = set(
        orbit(LEFT_DOUBLE_SINGLETON)
        + orbit(REVERSED_DOUBLE_SINGLETON)
        + orbit(DOUBLE_ENDPOINT_ONE_FOUR)
    )
    entries = set(degree_ten_completion_entries())
    if entries != expected_entries or len(entries) != 12:
        raise AssertionError(("degree-ten orbit inventory", entries))
    pre_degree_ten = set(pre_degree_ten_quintic_entries())
    remaining = set(remaining_quintic_entries())
    if not entries.issubset(pre_degree_ten):
        raise AssertionError("degree-ten entries outside live family")
    if entries.intersection(remaining):
        raise AssertionError("degree-ten partition overlap")
    if len(entries) + len(remaining) != len(pre_degree_ten):
        raise AssertionError("degree-ten partition size")
    if set(coefficient_map()) != entries:
        raise AssertionError("degree-ten coefficient map")
    for entry in entries:
        if np.isclose(
            inserted_coefficients()[entry],
            previous_inserted_coefficients()[entry],
        ):
            raise AssertionError(("degree-ten theorem not inserted", entry))

    q4 = q4_scalar_row_maxima()
    expected_q4 = (
        (0.42283950617284005, 2.759259259259264, 2.9212962962963007),
        (0.9722222222222242, 5.999999999999999, 5.999999999999999),
        (0.05246913580246914, 2.5000000000000013, 2.5000000000000013),
    )
    if not np.allclose(q4, expected_q4, rtol=1e-13, atol=1e-14):
        raise AssertionError(("degree-ten q4 exact rows", q4))
    residual_worst = double_singleton_residual_stress()

    committed = (
        ROOT / "artifacts" / "q64_degree_ten_completion_row_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 degree-ten completion-row artifact")

    print(
        "q64 degree-ten completion-row insertion passed: "
        f"entries={result.closed_entries},"
        f"coefficients={result.left_double_singleton_coefficient:.12g}/"
        f"{result.reversed_double_singleton_coefficient:.12g}/"
        f"{result.double_endpoint_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.degree_ten_inserted.total:.12g},"
        f"margin={result.degree_ten_inserted.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        f"residual_stress={residual_worst:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
