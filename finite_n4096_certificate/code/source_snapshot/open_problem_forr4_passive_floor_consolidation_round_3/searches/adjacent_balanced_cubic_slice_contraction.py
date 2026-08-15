#!/usr/bin/env python3
"""Chain-aware bound for the second unforced balanced high-sector orbit.

For profile ``(1,1,3,5)`` and split ``(0,0,1,4)``, rows are ``(x,F)``
and columns are ``(a,b,E,e)``.  The exact occurrence tensor is

    H(a,b) M_13(b, {x} union E) M_35({x} union E, F union {e}).

The column-only Walsh link contributes ``1/q``.  The internally split cubic
link has a Schur factor equal to the square root of its exact fixed-pair
squared slice.  Completing the adjacent moment as a unit-feature Gram symbol
leaves only the quintic 4|1 distinct-label Schur factor.  This proves a
coefficient below the provisional ``1/q`` target at ``N=1024``.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

from scipy.optimize import minimize_scalar

from attenuation_promise_concentration import (
    extended_euclidean_promise_concentration,
)
from leading_balanced_disjointness_contraction import (
    disjointness_schur_factor,
    leading_balanced_coefficient,
    leading_balanced_orbit_entries,
)
from occupation_compatible_sector_optimization import (
    certificate,
    endpoint_singleton_slice_energies,
)
from repaired_open_profile_budget import coarse_open_completion_coefficients


DIMENSION = 1024
THRESHOLD = 1 / 3
PROFILE = (1, 1, 3, 5)
SPLIT = (0, 0, 1, 4)


Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class AdjacentBalancedContraction:
    dimension: int
    order: int
    cubic_pair_slice_energy: float
    cubic_schur_factor: float
    quintic_disjointness_factor: float
    generic_two_mask_coefficient: float
    coefficient: float
    provisional_coefficient: float
    optimal_beta: float
    optimized_total: float
    threshold_slack: float


def adjacent_balanced_orbit_entries() -> tuple[ProfileSplit, ...]:
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


def adjacent_balanced_coefficient(
    dimension: int = DIMENSION,
) -> float:
    """Return the proved arbitrary-diagonal coefficient for the orbit."""

    order = int(round(sqrt(dimension)))
    if order * order != dimension:
        raise ValueError(("square dimension required", dimension))
    cubic_pair_slice = endpoint_singleton_slice_energies(order)[2]
    return (
        disjointness_schur_factor(dimension, 4)
        * sqrt(cubic_pair_slice)
        / order
    )


def adjacent_balanced_contraction(
    dimension: int = DIMENSION,
) -> AdjacentBalancedContraction:
    """Insert both proved balanced-orbit bounds into the coarse ledger."""

    if dimension != DIMENSION:
        raise ValueError(("ledger calibrated at N=1024", dimension))
    order = int(round(sqrt(dimension)))
    cubic_pair_slice = endpoint_singleton_slice_energies(order)[2]
    quintic_factor = disjointness_schur_factor(dimension, 4)
    generic_two_mask = (
        disjointness_schur_factor(dimension, 2)
        * quintic_factor
        / order
    )
    coefficient = adjacent_balanced_coefficient(dimension)
    coefficients = coarse_open_completion_coefficients()
    provisional = coefficients[(PROFILE, SPLIT)]
    for entry in leading_balanced_orbit_entries():
        coefficients[entry] = leading_balanced_coefficient(dimension)
    for entry in adjacent_balanced_orbit_entries():
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
    return AdjacentBalancedContraction(
        dimension=dimension,
        order=order,
        cubic_pair_slice_energy=cubic_pair_slice,
        cubic_schur_factor=sqrt(cubic_pair_slice),
        quintic_disjointness_factor=quintic_factor,
        generic_two_mask_coefficient=generic_two_mask,
        coefficient=coefficient,
        provisional_coefficient=provisional,
        optimal_beta=float(optimum.x),
        optimized_total=optimized,
        threshold_slack=THRESHOLD - optimized,
    )


def main() -> None:
    result = adjacent_balanced_contraction()
    print(
        "adjacent balanced cubic-slice contraction: "
        f"N={result.dimension},"
        f"q={result.order},"
        f"cubic_pair_slice={result.cubic_pair_slice_energy:.15g},"
        f"cubic_factor={result.cubic_schur_factor:.15g},"
        f"quintic_factor={result.quintic_disjointness_factor:.15g},"
        f"generic_two_mask={result.generic_two_mask_coefficient:.15g},"
        f"coefficient={result.coefficient:.15g},"
        f"provisional={result.provisional_coefficient:.15g},"
        f"optimal_beta={result.optimal_beta:.15g},"
        f"optimized_total={result.optimized_total:.15g},"
        f"threshold_slack={result.threshold_slack:.15g}"
    )


if __name__ == "__main__":
    main()
