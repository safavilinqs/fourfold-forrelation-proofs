#!/usr/bin/env python3
"""Contraction for the seventh balanced high-sector orbit.

For profile (3,1,1,5) and split (0,1,1,3), rows are indexed by

    (b, c, F),       |F| = 3,

and columns by

    (Q, G),          |Q| = 3, |G| = 2.

The exact occurrence tensor is

    M_31(Q,b) H(b,c) M_15(c,F union G).

The complete split-quintic row is used as a Schur feature.  After removing
it, duplicate compression and the constant modulus of H(b,c) contribute
1/sqrt(N).  Columns of M_31 with the same support xor are proportional
Walsh columns, and their squared cubic amplitudes have total mass at most
one.  Schatten Holder therefore bounds the remaining arbitrary-law weighted
trace norm by one.
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
from occupation_compatible_sector_optimization import certificate
from repaired_open_profile_budget import coarse_open_completion_coefficients
from separated_balanced_endpoint_slice_contraction import (
    separated_balanced_coefficient,
    separated_balanced_orbit_entries,
)


DIMENSION = 1024
THRESHOLD = 1 / 3
PROFILE = (3, 1, 1, 5)
SPLIT = (0, 1, 1, 3)
TARGET_GATE = 0.0484819899411186

Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class WholeCubicQuinticTripleContraction:
    dimension: int
    order: int
    quintic_three_slice_energy: float
    quintic_schur_factor: float
    compressed_cubic_hadamard_factor: float
    coefficient: float
    provisional_coefficient: float
    acceptance_gate: float
    optimal_beta: float
    optimized_total: float
    threshold_slack: float
    next_unresolved_entries: tuple[ProfileSplit, ...]
    next_unresolved_contribution: float
    next_admissible_coefficient: float


def whole_cubic_quintic_triple_orbit_entries() -> tuple[ProfileSplit, ...]:
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


def quintic_three_slice_energy(order: int) -> float:
    """Exact maximum squared M_15 tail through a fixed triple.

    The maximizing triple lies in one hidden column.  The three terms count
    quintic extensions of multiplicity types 5, 4+1, and 3+2 after the
    exceptional xor classes have been separated.
    """

    q = order
    if q < 4 or q & (q - 1):
        raise ValueError(("power-of-two order at least four required", q))
    w0 = 1 / q**2
    w1 = 1 / (q**2 * (q - 1) ** 2)
    return (
        comb(q - 3, 2) * w0
        + q * (q - 1) * (w0 + (q - 4) * w1)
        + (q - 1) * comb(q, 2) * w1
    )


def whole_cubic_quintic_triple_coefficient(
    dimension: int = DIMENSION,
) -> float:
    """Return the proved arbitrary-diagonal coefficient for the orbit."""

    order = int(round(sqrt(dimension)))
    if order * order != dimension:
        raise ValueError(("square dimension required", dimension))
    return sqrt(quintic_three_slice_energy(order) / dimension)


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


def whole_cubic_quintic_triple_contraction(
    dimension: int = DIMENSION,
) -> WholeCubicQuinticTripleContraction:
    """Insert the seventh theorem, rerank, and compute the next gate."""

    if dimension != DIMENSION:
        raise ValueError(("ledger calibrated at N=1024", dimension))
    q = int(round(sqrt(dimension)))
    slice_energy = quintic_three_slice_energy(q)
    coefficient = whole_cubic_quintic_triple_coefficient(dimension)

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
        (
            adjacent_balanced_row_orbit_entries(),
            adjacent_balanced_row_coefficient(),
        ),
        (
            whole_cubic_quintic_triple_orbit_entries(),
            coefficient,
        ),
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
    return WholeCubicQuinticTripleContraction(
        dimension=dimension,
        order=q,
        quintic_three_slice_energy=slice_energy,
        quintic_schur_factor=sqrt(slice_energy),
        compressed_cubic_hadamard_factor=1 / q,
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
    result = whole_cubic_quintic_triple_contraction()
    print(
        "whole-cubic quintic-triple contraction: "
        f"N={result.dimension},"
        f"q={result.order},"
        f"quintic_slice={result.quintic_three_slice_energy:.15g},"
        f"quintic_schur={result.quintic_schur_factor:.15g},"
        f"base_factor={result.compressed_cubic_hadamard_factor:.15g},"
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
