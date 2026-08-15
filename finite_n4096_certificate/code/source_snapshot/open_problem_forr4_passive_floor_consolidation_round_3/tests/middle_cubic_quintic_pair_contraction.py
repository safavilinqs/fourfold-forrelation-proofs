#!/usr/bin/env python3
"""Regression for the eighth balanced chained row-energy contraction."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys

import numpy as np
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from adjacent_cubic_quintic_mixed_orbit_q4 import (  # noqa: E402
    combined_link_moment,
)
from adjacent_cubic_quintic_orbit_witness import (  # noqa: E402
    parity_record_size,
    record_one_link_moment,
)
from middle_cubic_quintic_pair_contraction import (  # noqa: E402
    middle_cubic_quintic_pair_coefficient,
    middle_cubic_quintic_pair_contraction,
    middle_cubic_quintic_pair_orbit_entries,
    middle_cubic_quintic_pair_row_energy,
    middle_link_maxima,
)
from occupation_compatible_sector_optimization import (  # noqa: E402
    endpoint_quintic_singleton_slice_energies,
    endpoint_singleton_slice_energies,
)
from opposite_endpoint_orbit_scan import (  # noqa: E402
    cubic_weight,
    endpoint_moment,
    quintic_weight,
)


SEED = 2026071508


def transpose_support(
    order: int,
    support: tuple[int, ...],
) -> tuple[int, ...]:
    """Swap the physical row and column labels of a support."""

    return tuple(sorted((cell % order) * order + cell // order for cell in support))


def exact_q4_row_table() -> tuple[
    float,
    float,
    float,
    tuple[int, int],
    tuple[int, int],
    float,
    float,
]:
    """Compute every q=4 squared row energy by incidence products."""

    order = 4
    dimension = order * order
    cubics = tuple(combinations(range(dimension), 3))
    quintics = tuple(combinations(range(dimension), 5))
    pairs = tuple(combinations(range(dimension), 2))
    pair_index = {pair: index for index, pair in enumerate(pairs)}

    cubic_weights = np.asarray(
        [cubic_weight(transpose_support(order, cubic), order) for cubic in cubics]
    )
    quintic_weights = np.asarray(
        [quintic_weight(quintic, order) for quintic in quintics]
    )
    weighted = {
        1: np.zeros((len(cubics), len(quintics))),
        3: np.zeros((len(cubics), len(quintics))),
    }
    middle_maxima = {1: 0.0, 3: 0.0}
    for cubic_index, cubic in enumerate(cubics):
        if cubic_weights[cubic_index] == 0:
            continue
        record = parity_record_size(order, cubic, axis=1)
        for quintic_index, quintic in enumerate(quintics):
            if quintic_weights[quintic_index] == 0:
                continue
            moment = combined_link_moment(order, cubic, quintic)
            middle_maxima[record] = max(
                middle_maxima[record],
                abs(moment),
            )
            weighted[record][cubic_index, quintic_index] = (
                cubic_weights[cubic_index] ** 2
                * moment**2
                * quintic_weights[quintic_index] ** 2
                / dimension
            )

    cubic_incidence = sparse.lil_matrix((len(pairs), len(cubics)))
    for cubic_index, cubic in enumerate(cubics):
        for pair in combinations(cubic, 2):
            cubic_incidence[pair_index[pair], cubic_index] = 1
    quintic_incidence = sparse.lil_matrix((len(quintics), len(pairs)))
    for quintic_index, quintic in enumerate(quintics):
        for pair in combinations(quintic, 2):
            quintic_incidence[quintic_index, pair_index[pair]] = 1
    cubic_incidence = cubic_incidence.tocsr()
    quintic_incidence = quintic_incidence.tocsr()

    record_tables = {
        record: cubic_incidence @ table @ quintic_incidence
        for record, table in weighted.items()
    }
    total = record_tables[1] + record_tables[3]
    maximum_index = np.unravel_index(int(np.argmax(total)), total.shape)
    return (
        float(total[maximum_index]),
        float(record_tables[1].max()),
        float(record_tables[3].max()),
        pairs[maximum_index[0]],
        pairs[maximum_index[1]],
        middle_maxima[1],
        middle_maxima[3],
    )


def exact_q8_representative_row() -> tuple[float, float]:
    """Evaluate the q=8 horizontal/vertical row found by q=4 symmetry."""

    q = 8
    dimension = q * q
    fixed_cubic_pair = (0, 1)
    fixed_quintic_pair = (0, q)
    cubic_fixed = set(fixed_cubic_pair)
    quintic_fixed = set(fixed_quintic_pair)
    cubics = tuple(
        tuple(sorted(fixed_cubic_pair + (cell,)))
        for cell in range(dimension)
        if cell not in cubic_fixed
    )
    quintics = tuple(
        tuple(sorted(fixed_quintic_pair + added))
        for added in combinations(
            tuple(cell for cell in range(dimension) if cell not in quintic_fixed),
            3,
        )
    )
    energies = {1: 0.0, 3: 0.0}
    for cubic in cubics:
        cubic_factor = cubic_weight(transpose_support(q, cubic), q) ** 2
        if cubic_factor == 0:
            continue
        record = parity_record_size(q, cubic, axis=1)
        for quintic in quintics:
            quintic_factor = quintic_weight(quintic, q) ** 2
            if quintic_factor == 0:
                continue
            moment = combined_link_moment(q, cubic, quintic)
            energies[record] += cubic_factor * moment**2 * quintic_factor / dimension
    return energies[1], energies[3]


def exact_q4_endpoint_factors() -> tuple[float, float]:
    """Check the two separated endpoint sums used by the theorem."""

    q = 4
    dimension = q * q
    pairs = tuple(combinations(range(dimension), 2))
    maximum_cubic = 0.0
    for fixed in pairs:
        fixed_set = set(fixed)
        energy = 0.0
        for cell in range(dimension):
            if cell in fixed_set:
                continue
            cubic = tuple(sorted(fixed + (cell,)))
            for singleton in range(dimension):
                energy += (
                    record_one_link_moment(
                        q,
                        (singleton,),
                        cubic,
                    )
                    ** 2
                )
        maximum_cubic = max(maximum_cubic, energy)

    maximum_quintic = 0.0
    for fixed in pairs:
        fixed_set = set(fixed)
        available = tuple(cell for cell in range(dimension) if cell not in fixed_set)
        for singleton in range(dimension):
            energy = 0.0
            for added in combinations(available, 3):
                quintic = tuple(sorted(fixed + added))
                energy += (
                    endpoint_moment(
                        quintic,
                        singleton,
                        q,
                        5,
                        False,
                    )
                    ** 2
                )
            maximum_quintic = max(maximum_quintic, energy)
    return maximum_cubic, maximum_quintic


def direct_sparse_tensor_checks(q4_row_coefficient: float) -> float:
    """Stress exact q=4 target submatrices under correlated laws."""

    rng = np.random.default_rng(SEED)
    q = 4
    dimension = q * q
    pairs = tuple(combinations(range(dimension), 2))
    triples = tuple(combinations(range(dimension), 3))
    worst = 0.0
    for _ in range(12):
        rows = tuple(
            (
                pairs[int(rng.integers(len(pairs)))],
                pairs[int(rng.integers(len(pairs)))],
                int(rng.integers(dimension)),
            )
            for _ in range(26)
        )
        columns = tuple(
            (
                int(rng.integers(dimension)),
                int(rng.integers(dimension)),
                triples[int(rng.integers(len(triples)))],
            )
            for _ in range(40)
        )
        tensor = np.zeros((len(rows), len(columns)))
        for row_index, (fixed_cubic, fixed_quintic, final) in enumerate(rows):
            cubic_fixed = set(fixed_cubic)
            quintic_fixed = set(fixed_quintic)
            for column_index, (first, cell, triple) in enumerate(columns):
                if cell in cubic_fixed or quintic_fixed.intersection(triple):
                    continue
                cubic = tuple(sorted(fixed_cubic + (cell,)))
                quintic = tuple(sorted(fixed_quintic + triple))
                tensor[row_index, column_index] = (
                    record_one_link_moment(q, (first,), cubic)
                    * combined_link_moment(q, cubic, quintic)
                    * endpoint_moment(quintic, final, q, 5, False)
                )
        row_law = rng.dirichlet(np.ones(len(rows)))
        column_law = rng.dirichlet(np.ones(len(columns)))
        weighted = np.sqrt(row_law)[:, None] * tensor * np.sqrt(column_law)[None, :]
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > q4_row_coefficient * (1 + 5e-12):
            raise AssertionError(("sparse target tensor", nuclear, q4_row_coefficient))
        worst = max(worst, nuclear / q4_row_coefficient)
    return worst


def main() -> None:
    (
        q4_row,
        q4_record_one,
        q4_record_three,
        maximizing_cubic_pair,
        maximizing_quintic_pair,
        q4_middle_one,
        q4_middle_three,
    ) = exact_q4_row_table()
    if not np.isclose(q4_row, 0.1153978052126202, atol=3e-14):
        raise AssertionError(("q4 maximum row", q4_row))
    if not np.isclose(q4_record_one, 0.006172839506172835, atol=3e-14):
        raise AssertionError(("q4 record one", q4_record_one))
    if not np.isclose(q4_record_three, 0.11093964334705084, atol=3e-14):
        raise AssertionError(("q4 record three", q4_record_three))
    if maximizing_cubic_pair[0] // 4 != maximizing_cubic_pair[1] // 4:
        raise AssertionError(("q4 cubic pair geometry", maximizing_cubic_pair))
    if maximizing_quintic_pair[0] % 4 != maximizing_quintic_pair[1] % 4:
        raise AssertionError(("q4 quintic pair geometry", maximizing_quintic_pair))
    if not np.isclose(q4_middle_one, 1 / 4, atol=2e-14):
        raise AssertionError(("q4 record-one maximum", q4_middle_one))
    # The unrestricted record-three maximum is 3/4 at q=4, but the final
    # endpoint-record-one constraint removes those maximizers.
    if not np.isclose(q4_middle_three, 1 / 4, atol=2e-14):
        raise AssertionError(
            ("q4 endpoint-compatible record-three maximum", q4_middle_three)
        )

    q8_one, q8_three = exact_q8_representative_row()
    if not np.isclose(q8_one, 4.981124408669225e-06, atol=3e-16):
        raise AssertionError(("q8 record one", q8_one))
    if not np.isclose(q8_three, 0.0002844991599360749, atol=3e-15):
        raise AssertionError(("q8 record three", q8_three))

    q4_cubic, q4_quintic = exact_q4_endpoint_factors()
    expected_cubic = 16 * endpoint_singleton_slice_energies(4)[2]
    expected_quintic = endpoint_quintic_singleton_slice_energies(4)[2]
    if not np.isclose(q4_cubic, expected_cubic, atol=3e-14):
        raise AssertionError(("q4 cubic endpoint factor", q4_cubic))
    if not np.isclose(q4_quintic, expected_quintic, atol=3e-14):
        raise AssertionError(("q4 quintic endpoint factor", q4_quintic))

    for order in (4, 8, 16, 32):
        record_one, record_three, universal = middle_link_maxima(order)
        if universal != max(record_one, record_three):
            raise AssertionError(("middle maximum", order))
    if not middle_link_maxima(32)[0] > middle_link_maxima(32)[1]:
        raise AssertionError("record-one maximum should dominate at q=32")

    orbit = middle_cubic_quintic_pair_orbit_entries()
    if len(orbit) != 4:
        raise AssertionError(("eighth target orbit", orbit))
    coefficient = middle_cubic_quintic_pair_coefficient()
    if not np.isclose(coefficient, 0.02852815229229594, atol=3e-15):
        raise AssertionError(("q32 coefficient", coefficient))
    if not np.isclose(
        middle_cubic_quintic_pair_row_energy(),
        0.00081385547321243,
        atol=3e-16,
    ):
        raise AssertionError("q32 row energy")

    result = middle_cubic_quintic_pair_contraction()
    if not result.coefficient < result.provisional_coefficient:
        raise AssertionError(("provisional improvement", result))
    if not result.coefficient < result.acceptance_gate:
        raise AssertionError(("eighth-orbit gate", result))
    if not np.isclose(
        result.optimized_total,
        0.3328775891310207,
        atol=4e-11,
    ):
        raise AssertionError(("eight-theorem ledger", result))
    if not result.threshold_slack > 0.000455:
        raise AssertionError(("eight-theorem slack", result))
    if result.next_unresolved_entries[0] != (
        (1, 1, 3, 5),
        (0, 0, 3, 2),
    ):
        raise AssertionError(("reranked next orbit", result))
    if not np.isclose(
        result.next_admissible_coefficient,
        0.04543218921334891,
        atol=4e-14,
    ):
        raise AssertionError(("next coefficient gate", result))

    sparse_worst = direct_sparse_tensor_checks(np.sqrt(q4_row))
    print(
        "middle cubic-quintic pair contraction passed: "
        f"q4_row={q4_row:.12g},"
        f"q8_row={q8_one + q8_three:.12g},"
        f"q32_row={result.row_energy_bound:.12g},"
        f"coefficient={result.coefficient:.12g},"
        f"ledger_total={result.optimized_total:.12g},"
        f"threshold_slack={result.threshold_slack:.12g},"
        f"next_admissible={result.next_admissible_coefficient:.12g},"
        f"sparse_ratio={sparse_worst:.12g}"
    )


if __name__ == "__main__":
    main()
