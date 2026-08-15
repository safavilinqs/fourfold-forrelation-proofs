#!/usr/bin/env python3
"""Chain-aware bound for the first unforced balanced high-sector orbit.

For profile ``(3,1,1,5)`` and split ``(0,1,0,4)``, rows are ``(b,F)``
and columns are ``(Q,c,e)``.  On distinct labels the occurrence tensor is

    M_31(Q,b) H(b,c) M_51(F union {e},c).

The quintic occurrence mask is handled as a Schur multiplier.  After that
factor is removed, the first two links collapse exactly to a duplicated
Hadamard matrix and contribute ``1/sqrt(N)``.  The resulting arbitrary-law
coefficient is

    [1-k/N + sqrt(k(1-k/N)(1-1/N))] / sqrt(N),

with ``k=4``.  The script also inserts this proved coefficient into the
current coarse-completion ledger; every other open ``1/q`` charge remains a
route-selection target rather than a theorem.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

from scipy.optimize import minimize_scalar

from attenuation_promise_concentration import (
    extended_euclidean_promise_concentration,
)
from occupation_compatible_sector_optimization import certificate
from repaired_open_profile_budget import coarse_open_completion_coefficients


DIMENSION = 1024
THRESHOLD = 1 / 3
PROFILE = (3, 1, 1, 5)
SPLIT = (0, 1, 0, 4)


Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class LeadingBalancedContraction:
    dimension: int
    selected_quintic_labels: int
    disjointness_schur_factor: float
    coefficient: float
    provisional_coefficient: float
    optimal_beta: float
    optimized_total: float
    threshold_slack: float


def disjointness_schur_factor(dimension: int, selected: int) -> float:
    """Return an explicit gamma_2 bound for ``1[e not in F]``.

    For a selected ``k``-set ``F`` and a singleton ``e``, put

        x_F = 1_F - (k/N) 1,    y_e = 1_e - (1/N) 1.

    Then ``1[e not in F] = alpha - <x_F,y_e>`` with
    ``alpha=1-k/N``.  Optimally balancing the scalar and centered-vector
    coordinates gives the factor below.
    """

    if not 0 < selected < dimension:
        raise ValueError(("nontrivial selected set required", dimension, selected))
    alpha = 1 - selected / dimension
    return alpha + sqrt(
        selected * alpha * (1 - 1 / dimension)
    )


def leading_balanced_orbit_entries() -> tuple[ProfileSplit, ...]:
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


def leading_balanced_coefficient(dimension: int = DIMENSION) -> float:
    """Return the proved arbitrary-diagonal coefficient for the orbit."""

    order = int(round(sqrt(dimension)))
    if order * order != dimension:
        raise ValueError(("square dimension required", dimension))
    return disjointness_schur_factor(dimension, 4) / order


def leading_balanced_contraction(
    dimension: int = DIMENSION,
) -> LeadingBalancedContraction:
    """Insert the proved orbit bound into the coarse completion ledger."""

    if dimension != DIMENSION:
        raise ValueError(("ledger calibrated at N=1024", dimension))
    coefficient = leading_balanced_coefficient(dimension)
    coefficients = coarse_open_completion_coefficients()
    provisional = coefficients[(PROFILE, SPLIT)]
    for entry in leading_balanced_orbit_entries():
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
    return LeadingBalancedContraction(
        dimension=dimension,
        selected_quintic_labels=4,
        disjointness_schur_factor=disjointness_schur_factor(dimension, 4),
        coefficient=coefficient,
        provisional_coefficient=provisional,
        optimal_beta=float(optimum.x),
        optimized_total=optimized,
        threshold_slack=THRESHOLD - optimized,
    )


def main() -> None:
    result = leading_balanced_contraction()
    print(
        "leading balanced disjointness contraction: "
        f"N={result.dimension},"
        f"k={result.selected_quintic_labels},"
        f"schur_factor={result.disjointness_schur_factor:.15g},"
        f"coefficient={result.coefficient:.15g},"
        f"provisional={result.provisional_coefficient:.15g},"
        f"optimal_beta={result.optimal_beta:.15g},"
        f"optimized_total={result.optimized_total:.15g},"
        f"threshold_slack={result.threshold_slack:.15g}"
    )


if __name__ == "__main__":
    main()
