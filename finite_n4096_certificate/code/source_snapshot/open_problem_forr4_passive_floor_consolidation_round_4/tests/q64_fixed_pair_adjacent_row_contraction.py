#!/usr/bin/env python3
"""Regression for the q64 fixed-pair adjacent-row theorem."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import numpy as np
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_fixed_pair_adjacent_row_contraction import (  # noqa: E402
    artifact_text,
    cubic_completion_counts,
    diagnostic,
    fixed_pair_adjacent_coefficient,
    fixed_pair_adjacent_energy_parts,
    fixed_pair_adjacent_entries,
)
from adjacent_cubic_quintic_orbit_witness import (  # noqa: E402
    exact_link_moments,
    parity_record_size,
    record_one_link_moment,
)


def enumerated_completion_extrema(
    order: int,
) -> tuple[tuple[int, int, int], tuple[float, float, float]]:
    dimension = order * order
    maximum_counts = [0, 0, 0]
    maximum_squared = [0.0, 0.0, 0.0]
    for fixed_pair in combinations(range(dimension), 2):
        counts = [0, 0, 0]
        squared = [0.0, 0.0, 0.0]
        for extra in range(dimension):
            if extra in fixed_pair:
                continue
            cubic = tuple(sorted(fixed_pair + (extra,)))
            if parity_record_size(order, cubic, axis=0) != 1:
                continue
            column_record = parity_record_size(order, cubic, axis=1)
            if column_record == 1:
                kind = 0
            elif len({cell // order for cell in cubic}) == 1:
                kind = 1
            else:
                kind = 2
            counts[kind] += 1
            moment_squared = record_one_link_moment(
                order, (0,), cubic
            ) ** 2
            squared[kind] = max(squared[kind], moment_squared)
        maximum_counts = [
            max(old, new)
            for old, new in zip(maximum_counts, counts, strict=True)
        ]
        maximum_squared = [
            max(old, new)
            for old, new in zip(maximum_squared, squared, strict=True)
        ]
    return tuple(maximum_counts), tuple(maximum_squared)


def exact_q4_maximum_row_energy() -> float:
    order = 4
    dimension = order * order
    moments = exact_link_moments(order)
    pairs = tuple(combinations(range(dimension), 2))
    pair_index = {pair: index for index, pair in enumerate(pairs)}

    cubic_rows: list[int] = []
    cubic_columns: list[int] = []
    for cubic_index, cubic in enumerate(moments.supports_three):
        for pair in combinations(cubic, 2):
            cubic_rows.append(cubic_index)
            cubic_columns.append(pair_index[pair])
    cubic_incidence = sparse.csr_matrix(
        (
            np.ones(len(cubic_rows)),
            (cubic_rows, cubic_columns),
        ),
        shape=(len(moments.supports_three), len(pairs)),
    )

    quintic_rows: list[int] = []
    quintic_columns: list[int] = []
    for quintic_index, quintic in enumerate(moments.supports_five):
        for pair in combinations(quintic, 2):
            quintic_rows.append(quintic_index)
            quintic_columns.append(pair_index[pair])
    quintic_incidence = sparse.csr_matrix(
        (
            np.ones(len(quintic_rows)),
            (quintic_rows, quintic_columns),
        ),
        shape=(len(moments.supports_five), len(pairs)),
    )

    tails = np.square(moments.moment_35) @ quintic_incidence
    maximum = 0.0
    for singleton in range(dimension):
        energies = cubic_incidence.T @ (
            np.square(moments.moment_13[singleton])[:, None] * tails
        )
        maximum = max(maximum, float(energies.max()))
    return maximum


def main() -> None:
    result = diagnostic()
    observed_discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.closed_entries,
        result.record_one_cubic_completions,
        result.record_three_high_cubic_completions,
        result.record_three_low_cubic_completions,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
        result.remaining_quintic_proxy_entries,
    )
    expected_discrete = (
        64,
        4096,
        16_384,
        4,
        126,
        62,
        3906,
        260,
        264,
        624,
        164,
    )
    if observed_discrete != expected_discrete:
        raise AssertionError(("fixed-pair discrete result", observed_discrete))

    observed = (
        result.record_one_tail_bound,
        result.record_three_same_pair_tail_bound,
        result.record_three_distinct_pair_tail_bound,
        result.record_one_row_energy_bound,
        result.record_three_row_energy_bound,
        result.total_row_energy_bound,
        result.coefficient,
        result.previous_routing.total,
        result.fixed_pair_inserted.total,
        result.fixed_pair_inserted.beta,
        result.fixed_pair_inserted.perron_upper,
        result.fixed_pair_inserted.promise_loss,
        result.fixed_pair_inserted.margin_to_one_third,
        result.margin_improvement,
        result.remaining_quintic_local_proxy.total,
        result.remaining_quintic_local_proxy.margin_to_one_third,
        result.proxy_reserve_after_declared_allowance,
    )
    expected = (
        36.394989350747025,
        0.3002190828737629,
        0.012456403767218151,
        0.0002820792205384039,
        0.0046164640719675845,
        0.004898543292505988,
        0.06998959417303395,
        0.33056335386676317,
        0.3292482541342017,
        0.7460950063917603,
        0.31197416557950725,
        0.01727408855469446,
        0.0040850791991315916,
        0.0013150997325614422,
        0.33106696650853623,
        0.0022663668247970836,
        0.0012663668247970836,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("fixed-pair numeric result", observed))
    if not np.isclose(
        fixed_pair_adjacent_coefficient() ** 2,
        sum(fixed_pair_adjacent_energy_parts()),
        rtol=1e-13,
    ):
        raise AssertionError("fixed-pair coefficient identity")

    for order in (4, 8):
        counts, squared = enumerated_completion_extrema(order)
        expected_counts = cubic_completion_counts(order)
        expected_squared = (
            1 / (order * order * (order - 1) ** 2),
            1 / order**2,
            1 / (order * order * (order - 1) ** 2),
        )
        if counts != expected_counts or not np.allclose(
            squared, expected_squared, rtol=1e-13, atol=1e-15
        ):
            raise AssertionError(
                ("fixed-pair completion classification", order, counts, squared)
            )

    exact_q4 = exact_q4_maximum_row_energy()
    if not np.isclose(exact_q4, 0.46566358024691473, rtol=1e-12):
        raise AssertionError(("fixed-pair exact q4 row", exact_q4))
    if exact_q4 > sum(fixed_pair_adjacent_energy_parts(4)) + 1e-12:
        raise AssertionError("q4 row exceeds analytic majorant")

    if len(fixed_pair_adjacent_entries()) != 4:
        raise AssertionError("fixed-pair orbit size")
    committed = (
        ROOT / "artifacts" / "q64_fixed_pair_adjacent_row_contraction.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 fixed-pair adjacent artifact")

    print(
        "q64 fixed-pair adjacent row contraction passed: "
        f"entries={result.closed_entries},"
        f"coefficient={result.coefficient:.12g},"
        f"total={result.fixed_pair_inserted.total:.12g},"
        f"margin={result.fixed_pair_inserted.margin_to_one_third:.12g},"
        f"proxy_margin={result.remaining_quintic_local_proxy.margin_to_one_third:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
