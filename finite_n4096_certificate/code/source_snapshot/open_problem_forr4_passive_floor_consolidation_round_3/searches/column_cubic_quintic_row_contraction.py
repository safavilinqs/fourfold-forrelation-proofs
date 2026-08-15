#!/usr/bin/env python3
"""Row-energy contraction for the fifth balanced high-sector orbit.

For profile ``(3,1,1,5)`` and split ``(0,0,1,4)``, rows are ``(c,F)``
and columns are ``(Q,b,e)``.  The exact occurrence tensor is

    M_31(Q,b) H(b,c) M_15(c,F union {e}).

The quintic ``4|1`` factor is used as its own Schur feature.  Its exact
maximum row energy is ``1-4/N``.  After removing it, duplicate compression
leaves an ``N``-row matrix whose entries have modulus at most ``1/N``.
Rank--Frobenius therefore contributes ``1/sqrt(N)`` for arbitrary diagonal
laws.  The resulting coefficient is slightly below the provisional
``1/sqrt(N)`` charge at ``N=1024``.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

from scipy.optimize import brentq, minimize_scalar

from adjacent_balanced_cubic_slice_contraction import (
    adjacent_balanced_coefficient,
    adjacent_balanced_orbit_entries,
)
from attenuation_promise_concentration import (
    extended_euclidean_promise_concentration,
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
    disjointness_schur_factor,
    leading_balanced_coefficient,
    leading_balanced_orbit_entries,
)
from occupation_compatible_sector_optimization import (
    certificate,
    endpoint_quintic_singleton_slice_energies,
)
from repaired_open_profile_budget import coarse_open_completion_coefficients
from separated_balanced_endpoint_slice_contraction import (
    separated_balanced_coefficient,
    separated_balanced_orbit_entries,
)


DIMENSION = 1024
THRESHOLD = 1 / 3
PROFILE = (3, 1, 1, 5)
SPLIT = (0, 0, 1, 4)

Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class ColumnCubicQuinticContraction:
    dimension: int
    order: int
    quintic_four_slice_energy: float
    quintic_row_schur_factor: float
    base_rank_frobenius_factor: float
    generic_disjointness_coefficient: float
    coefficient: float
    provisional_coefficient: float
    acceptance_gate: float
    optimal_beta: float
    optimized_total: float
    threshold_slack: float
    next_unresolved_entries: tuple[ProfileSplit, ...]
    next_unresolved_contribution: float
    next_admissible_coefficient: float


def column_cubic_quintic_orbit_entries() -> tuple[ProfileSplit, ...]:
    """Return the complement/reversal orbit of the target split."""

    complement = tuple(
        degree - selected
        for degree, selected in zip(PROFILE, SPLIT, strict=True)
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


def column_cubic_quintic_coefficient(
    dimension: int = DIMENSION,
) -> float:
    """Return the proved arbitrary-diagonal coefficient for the orbit."""

    order = int(round(sqrt(dimension)))
    if order * order != dimension:
        raise ValueError(("square dimension required", dimension))
    quintic_row_energy = endpoint_quintic_singleton_slice_energies(order)[4]
    return sqrt(quintic_row_energy / dimension)


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


def column_cubic_quintic_contraction(
    dimension: int = DIMENSION,
) -> ColumnCubicQuinticContraction:
    """Insert the fifth theorem, rerank, and compute the next gate."""

    if dimension != DIMENSION:
        raise ValueError(("ledger calibrated at N=1024", dimension))
    order = int(round(sqrt(dimension)))
    quintic_row_energy = endpoint_quintic_singleton_slice_energies(order)[4]
    coefficient = column_cubic_quintic_coefficient(dimension)
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
        (column_cubic_quintic_orbit_entries(), coefficient),
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
    proved_orbits = tuple(
        frozenset(orbit) for orbit, _ in proved_orbit_values
    )
    unresolved = []
    for orbit in symmetry_orbits(dose_six_relevant_entries()):
        if frozenset(orbit) in proved_orbits:
            continue
        if not all(
            abs(base_coefficients[entry] - 1 / order) < 1e-14
            for entry in orbit
        ):
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
        1 / order,
        0.2,
        xtol=1e-14,
    )
    return ColumnCubicQuinticContraction(
        dimension=dimension,
        order=order,
        quintic_four_slice_energy=quintic_row_energy,
        quintic_row_schur_factor=sqrt(quintic_row_energy),
        base_rank_frobenius_factor=1 / order,
        generic_disjointness_coefficient=(
            disjointness_schur_factor(dimension, 4) / order
        ),
        coefficient=coefficient,
        provisional_coefficient=provisional,
        acceptance_gate=0.0542506297760259,
        optimal_beta=optimal_beta,
        optimized_total=optimized_total,
        threshold_slack=THRESHOLD - optimized_total,
        next_unresolved_entries=next_entries,
        next_unresolved_contribution=next_contribution,
        next_admissible_coefficient=next_admissible,
    )


def main() -> None:
    result = column_cubic_quintic_contraction()
    print(
        "column-cubic quintic-row contraction: "
        f"N={result.dimension},q={result.order},"
        f"quintic_row_energy={result.quintic_four_slice_energy:.15g},"
        f"quintic_schur={result.quintic_row_schur_factor:.15g},"
        f"base_factor={result.base_rank_frobenius_factor:.15g},"
        f"generic={result.generic_disjointness_coefficient:.15g},"
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
