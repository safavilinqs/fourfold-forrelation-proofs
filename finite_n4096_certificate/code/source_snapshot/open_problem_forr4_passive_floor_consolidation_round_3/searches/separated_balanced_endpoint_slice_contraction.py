#!/usr/bin/env python3
"""Exact-slice bound for the third balanced high-sector orbit.

For profile ``(3,1,1,5)`` and split ``(1,0,0,4)``, rows are ``(x,F)``
and columns are ``(E,b,c,e)``.  The exact occurrence tensor is

    M_31({x} union E,b) H(b,c) M_15(c,F union {e}).

The internally split cubic contributes the square root of its exact
fixed-pair slice.  Multiplying the central Walsh link into M_15 isolates the
endpoint quintic scalar v_5(F union {e}); the exact fixed-four-set moment
slice controls the resulting scalar/phase multiplier without the generic
4|1 distinct-label factor.  The residual Walsh phase is a unit-feature Schur
multiplier.  This gives a coefficient below the failing generic two-mask
bound at N=1024.
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
from leading_balanced_disjointness_contraction import (
    disjointness_schur_factor,
    leading_balanced_coefficient,
    leading_balanced_orbit_entries,
)
from occupation_compatible_sector_optimization import (
    certificate,
    endpoint_quintic_singleton_slice_energies,
    endpoint_singleton_slice_energies,
)
from repaired_open_profile_budget import coarse_open_completion_coefficients


DIMENSION = 1024
THRESHOLD = 1 / 3
PROFILE = (3, 1, 1, 5)
SPLIT = (1, 0, 0, 4)


Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class SeparatedBalancedContraction:
    dimension: int
    order: int
    cubic_pair_slice_energy: float
    cubic_schur_factor: float
    quintic_four_slice_energy: float
    collapsed_quintic_factor: float
    generic_two_mask_coefficient: float
    coefficient: float
    provisional_coefficient: float
    optimal_beta: float
    optimized_total: float
    threshold_slack: float
    next_unresolved_entries: tuple[ProfileSplit, ...]
    next_unresolved_contribution: float
    next_admissible_coefficient: float
    following_unresolved_entries: tuple[ProfileSplit, ...]
    following_unresolved_contribution: float


def separated_balanced_orbit_entries() -> tuple[ProfileSplit, ...]:
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


def separated_balanced_coefficient(
    dimension: int = DIMENSION,
) -> float:
    """Return the proved arbitrary-diagonal coefficient for the orbit."""

    order = int(round(sqrt(dimension)))
    if order * order != dimension:
        raise ValueError(("square dimension required", dimension))
    cubic_pair_slice = endpoint_singleton_slice_energies(order)[2]
    quintic_four_slice = endpoint_quintic_singleton_slice_energies(order)[4]
    return sqrt(cubic_pair_slice * quintic_four_slice)


def separated_balanced_contraction(
    dimension: int = DIMENSION,
) -> SeparatedBalancedContraction:
    """Insert all three proved balanced-orbit bounds into the coarse ledger."""

    if dimension != DIMENSION:
        raise ValueError(("ledger calibrated at N=1024", dimension))
    order = int(round(sqrt(dimension)))
    cubic_pair_slice = endpoint_singleton_slice_energies(order)[2]
    quintic_four_slice = endpoint_quintic_singleton_slice_energies(order)[4]
    generic_two_mask = (
        disjointness_schur_factor(dimension, 2)
        * disjointness_schur_factor(dimension, 4)
        / order
    )
    coefficient = separated_balanced_coefficient(dimension)
    base_coefficients = coarse_open_completion_coefficients()
    coefficients = dict(base_coefficients)
    provisional = base_coefficients[(PROFILE, SPLIT)]
    for entry in leading_balanced_orbit_entries():
        coefficients[entry] = leading_balanced_coefficient(dimension)
    for entry in adjacent_balanced_orbit_entries():
        coefficients[entry] = adjacent_balanced_coefficient(dimension)
    for entry in separated_balanced_orbit_entries():
        coefficients[entry] = coefficient

    def total(beta: float) -> float:
        ledger = certificate(
            beta=beta,
            profile_split_coefficients=coefficients,
        )
        promise = extended_euclidean_promise_concentration(
            dimension,
            beta,
        )
        return ledger.supporting_upper + promise.two_hypothesis_loss

    optimum = minimize_scalar(
        total,
        bounds=(0.75, 0.81),
        method="bounded",
        options={"xatol": 1e-13},
    )
    optimized = float(optimum.fun)
    ledger = certificate(
        beta=float(optimum.x),
        profile_split_coefficients=coefficients,
    )
    weights = dict(ledger.occupation_weights)
    proved_orbits = tuple(
        frozenset(orbit)
        for orbit in (
            leading_balanced_orbit_entries(),
            adjacent_balanced_orbit_entries(),
            separated_balanced_orbit_entries(),
        )
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
                float(optimum.x),
                weights,
            )
            * coefficients[(profile, split)]
            for profile, split in orbit
        )
        unresolved.append((contribution, orbit))
    unresolved.sort(reverse=True)
    next_contribution, next_entries = unresolved[0]
    following_contribution, following_entries = unresolved[1]

    def total_with_next_coefficient(value: float) -> float:
        trial_coefficients = dict(coefficients)
        for entry in next_entries:
            trial_coefficients[entry] = value

        def trial_total(beta: float) -> float:
            trial_ledger = certificate(
                beta=beta,
                profile_split_coefficients=trial_coefficients,
            )
            trial_promise = extended_euclidean_promise_concentration(
                dimension,
                beta,
            )
            return (
                trial_ledger.supporting_upper
                + trial_promise.two_hypothesis_loss
            )

        trial_optimum = minimize_scalar(
            trial_total,
            bounds=(0.75, 0.81),
            method="bounded",
            options={"xatol": 1e-13},
        )
        return float(trial_optimum.fun)

    next_admissible = brentq(
        lambda value: total_with_next_coefficient(value) - THRESHOLD,
        1 / order,
        0.2,
        xtol=1e-14,
    )
    return SeparatedBalancedContraction(
        dimension=dimension,
        order=order,
        cubic_pair_slice_energy=cubic_pair_slice,
        cubic_schur_factor=sqrt(cubic_pair_slice),
        quintic_four_slice_energy=quintic_four_slice,
        collapsed_quintic_factor=sqrt(quintic_four_slice),
        generic_two_mask_coefficient=generic_two_mask,
        coefficient=coefficient,
        provisional_coefficient=provisional,
        optimal_beta=float(optimum.x),
        optimized_total=optimized,
        threshold_slack=THRESHOLD - optimized,
        next_unresolved_entries=next_entries,
        next_unresolved_contribution=next_contribution,
        next_admissible_coefficient=next_admissible,
        following_unresolved_entries=following_entries,
        following_unresolved_contribution=following_contribution,
    )


def main() -> None:
    result = separated_balanced_contraction()
    print(
        "separated balanced endpoint-slice contraction: "
        f"N={result.dimension},"
        f"q={result.order},"
        f"cubic_pair_slice={result.cubic_pair_slice_energy:.15g},"
        f"cubic_factor={result.cubic_schur_factor:.15g},"
        f"quintic_four_slice={result.quintic_four_slice_energy:.15g},"
        f"collapsed_quintic_factor={result.collapsed_quintic_factor:.15g},"
        f"generic_two_mask={result.generic_two_mask_coefficient:.15g},"
        f"coefficient={result.coefficient:.15g},"
        f"provisional={result.provisional_coefficient:.15g},"
        f"optimal_beta={result.optimal_beta:.15g},"
        f"optimized_total={result.optimized_total:.15g},"
        f"threshold_slack={result.threshold_slack:.15g},"
        f"next_unresolved={result.next_unresolved_entries[0]},"
        f"next_contribution={result.next_unresolved_contribution:.15g},"
        f"next_admissible={result.next_admissible_coefficient:.15g},"
        f"next_admissible_over_provisional="
        f"{result.next_admissible_coefficient / result.provisional_coefficient:.15g},"
        f"following_unresolved={result.following_unresolved_entries[0]},"
        f"following_contribution="
        f"{result.following_unresolved_contribution:.15g}"
    )


if __name__ == "__main__":
    main()
