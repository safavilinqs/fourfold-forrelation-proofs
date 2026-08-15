#!/usr/bin/env python3
"""Row-slice contraction for the sixth balanced high-sector orbit.

For profile ``(1,1,3,5)`` and split ``(0,1,1,3)``, rows are indexed by

    (b, x, F),       |F| = 3,

and columns by

    (a, E, G),       |E| = |G| = 2.

With ``C={x} union E`` and ``S=F union G``, the exact occurrence kernel is

    H(a,b) M_13(b,C) M_35(C,S).

The intended arbitrary-law bound factors the complete ``M_13 M_35`` row as
a Schur feature.  The matrix left after this factor is the repeated Hadamard
matrix, whose weighted trace norm is at most the geometric mean of the row
and column masses.  Thus the target coefficient is the square root of the
largest chain row energy

    sum_{E,G} |M_13(b,{x} union E) M_35({x} union E,F union G)|^2.

The general-q proof separates the record-one and record-three sectors of
``M_35``.  Complete q=4 rows and selected q=8 rows protect the support
geometry, orthogonality identity, and normalization.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import comb, sqrt
import argparse

import numpy as np
from scipy import sparse
from scipy.optimize import brentq, minimize_scalar

from adjacent_cubic_quintic_mixed_orbit_q4 import combined_link_moment
from adjacent_cubic_quintic_orbit_witness import (
    parity_record_size,
    record_one_link_moment,
)
from adjacent_balanced_cubic_slice_contraction import (
    adjacent_balanced_coefficient,
    adjacent_balanced_orbit_entries,
)
from attenuation_promise_concentration import (
    extended_euclidean_promise_concentration,
)
from column_cubic_quintic_row_contraction import (
    column_cubic_quintic_coefficient,
    column_cubic_quintic_orbit_entries,
)
from high_degree_record_incidence_frontier import (
    dose_six_relevant_entries,
    split_perron_sensitivity,
    symmetry_orbits,
)
from internal_singleton_shared_law_contraction import (
    internal_singleton_coefficient,
    internal_singleton_orbit_entries,
)
from leading_balanced_disjointness_contraction import (
    leading_balanced_coefficient,
    leading_balanced_orbit_entries,
)
from occupation_compatible_sector_optimization import certificate
from repaired_open_profile_budget import coarse_open_completion_coefficients
from separated_balanced_endpoint_slice_contraction import (
    separated_balanced_coefficient,
    separated_balanced_orbit_entries,
)


PROFILE = (1, 1, 3, 5)
SPLIT = (0, 1, 1, 3)
TARGET_ORDER = 32
TARGET_GATE = 0.0570749885142
DIMENSION = TARGET_ORDER**2
THRESHOLD = 1 / 3

Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class ExactRowSlice:
    order: int
    fixed_cubic_cell: int
    rows: int
    cubic_extensions: int
    quintic_fixed_triples: int
    maximum_row_energy: float
    maximum_record_one_energy: float
    maximum_record_three_energy: float
    maximum_coefficient: float
    maximizing_singleton: int
    maximizing_triple: tuple[int, int, int]
    maximizing_record_one: float
    maximizing_record_three: float


@dataclass(frozen=True)
class FixedRowEnergy:
    order: int
    singleton: int
    fixed_cubic_cell: int
    fixed_quintic_triple: tuple[int, int, int]
    cubic_extensions: int
    quintic_extensions: int
    record_one_energy: float
    record_three_energy: float
    total_energy: float
    coefficient: float


@dataclass(frozen=True)
class AdjacentBalancedRowContraction:
    dimension: int
    order: int
    record_one_row_energy_bound: float
    record_three_cubic_weight: float
    horizontal_record_three_tail: float
    two_one_record_three_tail_bound: float
    distinct_record_three_tail_bound: float
    record_three_row_energy_bound: float
    total_row_energy_bound: float
    coefficient: float
    provisional_coefficient: float
    acceptance_gate: float
    optimal_beta: float
    optimized_total: float
    threshold_slack: float
    next_unresolved_entries: tuple[ProfileSplit, ...]
    next_unresolved_contribution: float
    next_admissible_coefficient: float


def target_orbit_entries() -> tuple[ProfileSplit, ...]:
    """Return complement/reversal orbit of the target cut."""

    complement = tuple(
        degree - selected for degree, selected in zip(PROFILE, SPLIT, strict=True)
    )
    reverse = tuple(reversed(PROFILE))
    return tuple(
        sorted(
            {
                (PROFILE, SPLIT),
                (PROFILE, complement),
                (reverse, tuple(reversed(SPLIT))),
                (reverse, tuple(reversed(complement))),
            }
        )
    )


def record_one_row_energy_bound(order: int) -> float:
    """Universal record-one part of the complete chain row energy.

    A cubic with record one on both link axes is one of the
    ``3(q-1)^2`` L-shapes through the fixed cell.  Its preceding singleton
    moment has squared modulus ``1/[q^2(q-1)^2]``.  A fixed three-cell subset
    has at most ``(q-2)(3q^2+q-6)/2`` record-one quintic extensions, and the
    cubic--quintic entry is at most
    ``(q+2)/[q(q-1)(q-2)]``.
    """

    q = order
    if q < 4 or q & (q - 1):
        raise ValueError(("power-of-two order at least four required", q))
    quintic_extensions = (q - 2) * (3 * q * q + q - 6) / 2
    return 3 * quintic_extensions * (q + 2) ** 2 / (q**4 * (q - 1) ** 2 * (q - 2) ** 2)


def record_three_cubic_weight(order: int) -> float:
    """Sum of ``M_13^2`` over record-three cubics containing one cell."""

    q = order
    same_row = comb(q - 1, 2)
    total = (q - 1) * (q - 2) * (3 * q - 2) / 2
    return same_row / q**2 + (total - same_row) / (q**2 * (q - 1) ** 2)


def record_three_tail_bounds(order: int) -> tuple[float, float, float]:
    """M_35 output bounds for the three row patterns of the fixed triple.

    The return order is ``(3)``, ``(2+1)``, and ``(1+1+1)``.  For a
    three-in-one-row triple, Walsh orthogonality makes the exact squared tail
    ``3/[(q-1)(q-2)]``.  The other two cases use exact extension counts and
    the no-even/one-even record-three entry bounds.
    """

    q = order
    denominator = comb(q, 3) ** 2
    horizontal = 3 / ((q - 1) * (q - 2))

    two_one_no_even = q * (q - 2) ** 2
    two_one_one_even = q**2 * comb(q - 2, 2)
    two_one = two_one_no_even / denominator + 9 * two_one_one_even / (
        (q - 3) ** 2 * denominator
    )

    distinct_no_even = 3 * comb(q - 1, 2)
    distinct_one_even = (q - 3) * comb(q, 2) + 3 * (q - 1) * (q - 3) * q
    distinct = distinct_no_even / denominator + 9 * distinct_one_even / (
        (q - 3) ** 2 * denominator
    )
    return horizontal, two_one, distinct


def adjacent_balanced_row_coefficient(
    dimension: int = DIMENSION,
) -> float:
    """Return the proved arbitrary-diagonal coefficient for the target."""

    order = int(round(sqrt(dimension)))
    if order * order != dimension:
        raise ValueError(("square dimension required", dimension))
    record_one = record_one_row_energy_bound(order)
    record_three = record_three_cubic_weight(order) * max(
        record_three_tail_bounds(order)
    )
    return sqrt(record_one + record_three)


def _optimize_total(
    coefficients: dict[ProfileSplit, float],
    dimension: int,
) -> tuple[float, float]:
    def total(beta: float) -> float:
        ledger = certificate(
            beta=beta,
            profile_split_coefficients=coefficients,
        )
        promise = extended_euclidean_promise_concentration(dimension, beta)
        return ledger.supporting_upper + promise.two_hypothesis_loss

    optimum = minimize_scalar(
        total,
        bounds=(0.75, 0.81),
        method="bounded",
        options={"xatol": 1e-13},
    )
    return float(optimum.x), float(optimum.fun)


def adjacent_balanced_row_contraction(
    dimension: int = DIMENSION,
) -> AdjacentBalancedRowContraction:
    """Insert the sixth theorem, rerank, and calculate the next gate."""

    if dimension != DIMENSION:
        raise ValueError(("ledger calibrated at N=1024", dimension))
    q = int(round(sqrt(dimension)))
    record_one = record_one_row_energy_bound(q)
    cubic_weight = record_three_cubic_weight(q)
    tails = record_three_tail_bounds(q)
    record_three = cubic_weight * max(tails)
    total_energy = record_one + record_three
    coefficient = sqrt(total_energy)

    base_coefficients = coarse_open_completion_coefficients()
    coefficients = dict(base_coefficients)
    provisional = base_coefficients[(PROFILE, SPLIT)]
    proved_orbit_values = (
        (leading_balanced_orbit_entries(), leading_balanced_coefficient()),
        (adjacent_balanced_orbit_entries(), adjacent_balanced_coefficient()),
        (
            separated_balanced_orbit_entries(),
            separated_balanced_coefficient(),
        ),
        (
            internal_singleton_orbit_entries(),
            internal_singleton_coefficient(),
        ),
        (
            column_cubic_quintic_orbit_entries(),
            column_cubic_quintic_coefficient(),
        ),
        (target_orbit_entries(), coefficient),
    )
    for orbit, value in proved_orbit_values:
        for entry in orbit:
            coefficients[entry] = value

    optimal_beta, optimized_total = _optimize_total(coefficients, dimension)
    ledger = certificate(
        beta=optimal_beta,
        profile_split_coefficients=coefficients,
    )
    weights = dict(ledger.occupation_weights)
    proved_orbits = tuple(frozenset(orbit) for orbit, _ in proved_orbit_values)
    unresolved = []
    for orbit in symmetry_orbits(dose_six_relevant_entries()):
        if frozenset(orbit) in proved_orbits:
            continue
        if not all(abs(base_coefficients[entry] - 1 / q) < 1e-14 for entry in orbit):
            continue
        contribution = sum(
            split_perron_sensitivity(
                profile,
                split,
                optimal_beta,
                weights,
            )
            * coefficients[(profile, split)]
            for profile, split in orbit
        )
        unresolved.append((contribution, orbit))
    unresolved.sort(reverse=True)
    next_contribution, next_entries = unresolved[0]

    def total_with_next(value: float) -> float:
        trial = dict(coefficients)
        for entry in next_entries:
            trial[entry] = value
        return _optimize_total(trial, dimension)[1]

    next_admissible = brentq(
        lambda value: total_with_next(value) - THRESHOLD,
        1 / q,
        0.3,
        xtol=1e-14,
    )
    return AdjacentBalancedRowContraction(
        dimension=dimension,
        order=q,
        record_one_row_energy_bound=record_one,
        record_three_cubic_weight=cubic_weight,
        horizontal_record_three_tail=tails[0],
        two_one_record_three_tail_bound=tails[1],
        distinct_record_three_tail_bound=tails[2],
        record_three_row_energy_bound=record_three,
        total_row_energy_bound=total_energy,
        coefficient=coefficient,
        provisional_coefficient=provisional,
        acceptance_gate=TARGET_GATE,
        optimal_beta=optimal_beta,
        optimized_total=optimized_total,
        threshold_slack=THRESHOLD - optimized_total,
        next_unresolved_entries=next_entries,
        next_unresolved_contribution=next_contribution,
        next_admissible_coefficient=next_admissible,
    )


def exact_q4_row_slice(
    *,
    fixed_cubic_cell: int = 0,
) -> ExactRowSlice:
    """Compute every q=4 row energy after translating ``x`` to one cell."""

    order = 4
    dimension = order**2
    if not 0 <= fixed_cubic_cell < dimension:
        raise ValueError(fixed_cubic_cell)
    available = tuple(cell for cell in range(dimension) if cell != fixed_cubic_cell)
    cubics = tuple(
        tuple(sorted((fixed_cubic_cell,) + pair)) for pair in combinations(available, 2)
    )
    triples = tuple(combinations(range(dimension), 3))
    quintics = tuple(combinations(range(dimension), 5))
    triple_index = {triple: index for index, triple in enumerate(triples)}

    middle_squared = np.asarray(
        [
            [
                record_one_link_moment(order, (singleton,), cubic) ** 2
                for cubic in cubics
            ]
            for singleton in range(dimension)
        ],
        dtype=float,
    )

    record_one = np.zeros((len(cubics), len(quintics)), dtype=float)
    record_three = np.zeros_like(record_one)
    for cubic_index, cubic in enumerate(cubics):
        record = parity_record_size(order, cubic, axis=1)
        target = record_one if record == 1 else record_three
        for quintic_index, quintic in enumerate(quintics):
            moment = combined_link_moment(order, cubic, quintic)
            target[cubic_index, quintic_index] = moment * moment

    incidence_rows: list[int] = []
    incidence_columns: list[int] = []
    for quintic_index, quintic in enumerate(quintics):
        for triple in combinations(quintic, 3):
            incidence_rows.append(quintic_index)
            incidence_columns.append(triple_index[triple])
    incidence = sparse.csr_matrix(
        (
            np.ones(len(incidence_rows), dtype=float),
            (incidence_rows, incidence_columns),
        ),
        shape=(len(quintics), len(triples)),
    )
    record_one_tail = record_one @ incidence
    record_three_tail = record_three @ incidence
    record_one_energy = middle_squared @ record_one_tail
    record_three_energy = middle_squared @ record_three_tail
    total = record_one_energy + record_three_energy

    maximum_index = np.unravel_index(int(np.argmax(total)), total.shape)
    singleton_index, triple_position = map(int, maximum_index)
    maximum_record_one = float(record_one_energy.max())
    maximum_record_three = float(record_three_energy.max())
    maximum = float(total[maximum_index])
    return ExactRowSlice(
        order=order,
        fixed_cubic_cell=fixed_cubic_cell,
        rows=dimension * len(triples),
        cubic_extensions=len(cubics),
        quintic_fixed_triples=len(triples),
        maximum_row_energy=maximum,
        maximum_record_one_energy=maximum_record_one,
        maximum_record_three_energy=maximum_record_three,
        maximum_coefficient=sqrt(maximum),
        maximizing_singleton=singleton_index,
        maximizing_triple=triples[triple_position],
        maximizing_record_one=float(record_one_energy[maximum_index]),
        maximizing_record_three=float(record_three_energy[maximum_index]),
    )


def exact_fixed_row_energy(
    order: int,
    singleton: int,
    fixed_cubic_cell: int,
    fixed_quintic_triple: tuple[int, int, int],
) -> FixedRowEnergy:
    """Direct exact row energy, intended for selected moderate-q rows."""

    dimension = order**2
    if len(set(fixed_quintic_triple)) != 3:
        raise ValueError(fixed_quintic_triple)
    if not all(0 <= cell < dimension for cell in fixed_quintic_triple):
        raise ValueError(fixed_quintic_triple)
    cubic_available = tuple(
        cell for cell in range(dimension) if cell != fixed_cubic_cell
    )
    cubics = tuple(
        tuple(sorted((fixed_cubic_cell,) + pair))
        for pair in combinations(cubic_available, 2)
    )
    fixed_set = set(fixed_quintic_triple)
    quintic_available = tuple(
        cell for cell in range(dimension) if cell not in fixed_set
    )
    quintics = tuple(
        tuple(sorted(fixed_quintic_triple + pair))
        for pair in combinations(quintic_available, 2)
    )
    record_energies = {1: 0.0, 3: 0.0}
    for cubic in cubics:
        middle = record_one_link_moment(order, (singleton,), cubic)
        if middle == 0:
            continue
        record = parity_record_size(order, cubic, axis=1)
        tail = sum(
            combined_link_moment(order, cubic, quintic) ** 2 for quintic in quintics
        )
        record_energies[record] += middle * middle * tail
    total = record_energies[1] + record_energies[3]
    return FixedRowEnergy(
        order=order,
        singleton=singleton,
        fixed_cubic_cell=fixed_cubic_cell,
        fixed_quintic_triple=fixed_quintic_triple,
        cubic_extensions=len(cubics),
        quintic_extensions=len(quintics),
        record_one_energy=record_energies[1],
        record_three_energy=record_energies[3],
        total_energy=total,
        coefficient=sqrt(total),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixed-order", type=int)
    parser.add_argument("--exact-q4", action="store_true")
    arguments = parser.parse_args()
    if arguments.fixed_order:
        order = arguments.fixed_order
        fixed = (order, 2 * order, 3 * order)
        result = exact_fixed_row_energy(order, 0, 0, fixed)
        print(
            "adjacent balanced fixed vertical row: "
            f"q={result.order},"
            f"b={result.singleton},"
            f"x={result.fixed_cubic_cell},"
            f"F={result.fixed_quintic_triple},"
            f"cubics={result.cubic_extensions},"
            f"quintics={result.quintic_extensions},"
            f"record_one={result.record_one_energy:.15g},"
            f"record_three={result.record_three_energy:.15g},"
            f"energy={result.total_energy:.15g},"
            f"coefficient={result.coefficient:.15g}"
        )
        return
    if not arguments.exact_q4:
        result = adjacent_balanced_row_contraction()
        print(
            "adjacent balanced row-slice contraction: "
            f"N={result.dimension},"
            f"q={result.order},"
            f"record_one={result.record_one_row_energy_bound:.15g},"
            f"record_three_weight={result.record_three_cubic_weight:.15g},"
            f"horizontal_tail={result.horizontal_record_three_tail:.15g},"
            f"two_one_tail={result.two_one_record_three_tail_bound:.15g},"
            f"distinct_tail={result.distinct_record_three_tail_bound:.15g},"
            f"record_three={result.record_three_row_energy_bound:.15g},"
            f"row_energy={result.total_row_energy_bound:.15g},"
            f"coefficient={result.coefficient:.15g},"
            f"provisional={result.provisional_coefficient:.15g},"
            f"gate={result.acceptance_gate:.15g},"
            f"optimal_beta={result.optimal_beta:.15g},"
            f"optimized_total={result.optimized_total:.15g},"
            f"threshold_slack={result.threshold_slack:.15g},"
            f"next_unresolved={result.next_unresolved_entries[0]},"
            f"next_contribution={result.next_unresolved_contribution:.15g},"
            f"next_admissible={result.next_admissible_coefficient:.15g}"
        )
        return
    result = exact_q4_row_slice()
    print(
        "adjacent balanced exact row slice: "
        f"q={result.order},"
        f"x={result.fixed_cubic_cell},"
        f"rows={result.rows},"
        f"cubics={result.cubic_extensions},"
        f"triples={result.quintic_fixed_triples},"
        f"max_energy={result.maximum_row_energy:.15g},"
        f"max_coefficient={result.maximum_coefficient:.15g},"
        f"max_record_one={result.maximum_record_one_energy:.15g},"
        f"max_record_three={result.maximum_record_three_energy:.15g},"
        f"argmax_b={result.maximizing_singleton},"
        f"argmax_F={result.maximizing_triple},"
        f"argmax_record_one={result.maximizing_record_one:.15g},"
        f"argmax_record_three={result.maximizing_record_three:.15g},"
        f"target_q={TARGET_ORDER},"
        f"target_gate={TARGET_GATE:.15g},"
        f"orbit={target_orbit_entries()}"
    )


if __name__ == "__main__":
    main()
