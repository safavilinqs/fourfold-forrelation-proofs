#!/usr/bin/env python3
"""Whole-cubic/middle-pair contraction for the ninth balanced orbit.

For profile (1,1,3,5) and split (0,0,3,2), rows are indexed by

    (C, F),       |C| = 3, |F| = 2,

and columns by

    (a, b, G),    |G| = 3.

With S=F union G, the exact occurrence tensor is

    H(a,b) M_13(b,C) M_35(C,S).

Extract a normalized complete M_35 row as a Schur feature.  The remaining
H M_13 chain collapses first by cubic XOR and then by the repeated a-label,
giving an additional exact 1/q.  Explicit fixed-pair extension counts bound
the record-one and record-three M_35 rows.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb, sqrt

from scipy.optimize import brentq, minimize_scalar

from adjacent_balanced_cubic_slice_contraction import (
    adjacent_balanced_coefficient,
    adjacent_balanced_orbit_entries,
)
from adjacent_balanced_row_slice_contraction import (
    adjacent_balanced_row_coefficient,
    target_orbit_entries as adjacent_balanced_row_orbit_entries,
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
from middle_cubic_quintic_pair_contraction import (
    middle_cubic_quintic_pair_coefficient,
    middle_cubic_quintic_pair_orbit_entries,
)
from occupation_compatible_sector_optimization import certificate
from repaired_open_profile_budget import coarse_open_completion_coefficients
from separated_balanced_endpoint_slice_contraction import (
    separated_balanced_coefficient,
    separated_balanced_orbit_entries,
)
from whole_cubic_quintic_triple_contraction import (
    whole_cubic_quintic_triple_coefficient,
    whole_cubic_quintic_triple_orbit_entries,
)


DIMENSION = 1024
THRESHOLD = 1 / 3
PROFILE = (1, 1, 3, 5)
SPLIT = (0, 0, 3, 2)
TARGET_GATE = 0.04543218921334891

Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class RecordThreePairCounts:
    same_row_no_even: int
    same_row_one_even: int
    distinct_rows_no_even: int
    distinct_rows_one_even: int


@dataclass(frozen=True)
class WholeCubicMiddlePairContraction:
    dimension: int
    order: int
    record_one_pair_incidence: int
    record_one_middle_maximum: float
    record_one_slice_bound: float
    record_one_coefficient: float
    record_three_counts: RecordThreePairCounts
    record_three_same_row_slice_bound: float
    record_three_distinct_rows_slice_bound: float
    record_three_coefficient: float
    coefficient: float
    provisional_coefficient: float
    acceptance_gate: float
    optimal_beta: float
    optimized_total: float
    threshold_slack: float
    next_unresolved_entries: tuple[ProfileSplit, ...]
    next_unresolved_contribution: float
    next_admissible_coefficient: float


def whole_cubic_middle_pair_orbit_entries() -> tuple[ProfileSplit, ...]:
    """Return the complement/reversal orbit of the target cut."""

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


def record_one_pair_incidence(order: int) -> int:
    """Maximum record-one quintic extensions through a fixed pair."""

    q = order
    if q < 4:
        raise ValueError(("order at least four required", q))
    # The maximum has both fixed cells in one row.  The three new cells
    # either lie in one row or have multiplicities 2+1.  Exactly one final
    # row must have odd multiplicity.
    return (
        comb(q - 2, 3)
        + (q - 1) * comb(q, 3)
        + (q - 2) * (q - 1) * comb(q, 2)
        + (q - 1) * q * (comb(q - 2, 2) + (q - 2) * comb(q, 2))
    )


def record_three_pair_counts(order: int) -> RecordThreePairCounts:
    """Count fixed-pair quintics with three odd row labels.

    ``no_even`` denotes row multiplicity 3+1+1.  ``one_even`` denotes
    2+1+1+1 and receives the extra injective-average suppression.
    """

    q = order
    if q < 4:
        raise ValueError(("order at least four required", q))
    same_no_even = comb(q - 1, 2) * (q - 2) * q**2
    same_one_even = comb(q - 1, 3) * q**3
    distinct_no_even = (q - 2) * comb(q, 3) + 2 * (q - 2) * q * comb(q - 1, 2)
    distinct_one_even = (q - 2) * q * (q - 3) * comb(q, 2) + 2 * (q - 1) * comb(
        q - 2, 2
    ) * q**2
    return RecordThreePairCounts(
        same_row_no_even=same_no_even,
        same_row_one_even=same_one_even,
        distinct_rows_no_even=distinct_no_even,
        distinct_rows_one_even=distinct_one_even,
    )


def record_sector_bounds(
    order: int,
) -> tuple[int, float, float, float, RecordThreePairCounts, float, float, float]:
    """Return the two M_35 slice bounds and final sector coefficients."""

    q = order
    incidence = record_one_pair_incidence(q)
    record_one_maximum = (q + 2) / (q * (q - 1) * (q - 2))
    record_one_slice = incidence * record_one_maximum**2
    # A record-one cubic compatible with M_13 is an L-shape, whose endpoint
    # amplitude is 1/(q-1).  The residual Walsh chain contributes 1/q.
    record_one_coefficient = sqrt(record_one_slice) / (q * (q - 1))

    counts = record_three_pair_counts(q)
    record_three_no_even_maximum = 1 / comb(q, 3)
    record_three_one_even_maximum = 3 / ((q - 3) * comb(q, 3))
    same_row_slice = (
        counts.same_row_no_even * record_three_no_even_maximum**2
        + counts.same_row_one_even * record_three_one_even_maximum**2
    )
    distinct_rows_slice = (
        counts.distinct_rows_no_even * record_three_no_even_maximum**2
        + counts.distinct_rows_one_even * record_three_one_even_maximum**2
    )
    # A three-record cubic compatible with M_13 has endpoint amplitude one.
    record_three_coefficient = sqrt(max(same_row_slice, distinct_rows_slice)) / q
    return (
        incidence,
        record_one_maximum,
        record_one_slice,
        record_one_coefficient,
        counts,
        same_row_slice,
        distinct_rows_slice,
        record_three_coefficient,
    )


def whole_cubic_middle_pair_coefficient(
    dimension: int = DIMENSION,
) -> float:
    """Return the proved arbitrary-diagonal coefficient for the orbit."""

    q = int(round(sqrt(dimension)))
    if q * q != dimension:
        raise ValueError(("square dimension required", dimension))
    bounds = record_sector_bounds(q)
    return max(bounds[3], bounds[7])


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


def whole_cubic_middle_pair_contraction(
    dimension: int = DIMENSION,
) -> WholeCubicMiddlePairContraction:
    """Insert the ninth theorem, rerank, and compute the next gate."""

    if dimension != DIMENSION:
        raise ValueError(("ledger calibrated at N=1024", dimension))
    q = int(round(sqrt(dimension)))
    (
        incidence,
        record_one_maximum,
        record_one_slice,
        record_one_coefficient,
        counts,
        same_row_slice,
        distinct_rows_slice,
        record_three_coefficient,
    ) = record_sector_bounds(q)
    coefficient = max(record_one_coefficient, record_three_coefficient)

    base_coefficients = coarse_open_completion_coefficients()
    coefficients = dict(base_coefficients)
    provisional = base_coefficients[(PROFILE, SPLIT)]
    proved_orbit_values = (
        (leading_balanced_orbit_entries(), leading_balanced_coefficient()),
        (adjacent_balanced_orbit_entries(), adjacent_balanced_coefficient()),
        (separated_balanced_orbit_entries(), separated_balanced_coefficient()),
        (internal_singleton_orbit_entries(), internal_singleton_coefficient()),
        (column_cubic_quintic_orbit_entries(), column_cubic_quintic_coefficient()),
        (
            adjacent_balanced_row_orbit_entries(),
            adjacent_balanced_row_coefficient(),
        ),
        (
            whole_cubic_quintic_triple_orbit_entries(),
            whole_cubic_quintic_triple_coefficient(),
        ),
        (
            middle_cubic_quintic_pair_orbit_entries(),
            middle_cubic_quintic_pair_coefficient(),
        ),
        (whole_cubic_middle_pair_orbit_entries(), coefficient),
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
    return WholeCubicMiddlePairContraction(
        dimension=dimension,
        order=q,
        record_one_pair_incidence=incidence,
        record_one_middle_maximum=record_one_maximum,
        record_one_slice_bound=record_one_slice,
        record_one_coefficient=record_one_coefficient,
        record_three_counts=counts,
        record_three_same_row_slice_bound=same_row_slice,
        record_three_distinct_rows_slice_bound=distinct_rows_slice,
        record_three_coefficient=record_three_coefficient,
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


def main() -> None:
    result = whole_cubic_middle_pair_contraction()
    print(
        "whole-cubic middle-pair contraction: "
        f"N={result.dimension},"
        f"q={result.order},"
        f"record_one_incidence={result.record_one_pair_incidence},"
        f"record_one_slice={result.record_one_slice_bound:.15g},"
        f"record_one_coefficient={result.record_one_coefficient:.15g},"
        f"record_three_same={result.record_three_same_row_slice_bound:.15g},"
        f"record_three_distinct={result.record_three_distinct_rows_slice_bound:.15g},"
        f"record_three_coefficient={result.record_three_coefficient:.15g},"
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


if __name__ == "__main__":
    main()
