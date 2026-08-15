#!/usr/bin/env python3
"""Double-endpoint row contraction for the tenth balanced orbit.

For profile (1,3,5,1) and split (0,2,3,0), rows are indexed by

    (E, F),       |E| = 2, |F| = 3,

and columns by

    (a, x, G, d), |x| = 1, |G| = 2.

With C=E union {x} and S=F union G, the exact occurrence tensor is

    M_13(a,C) M_35(C,S) M_51(S,d).

Extract the scalar completion row v_3(C) M_35(C,S) v_5(S) as a Schur
feature.  The two endpoint Walsh factors form a repeated, column-twisted
H_N tensor H_N base, whose arbitrary-diagonal weighted trace norm has
coefficient one.  Record-one incidence and record-three endpoint slices
bound the scalar row.
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
from occupation_compatible_sector_optimization import (
    certificate,
    endpoint_quintic_singleton_slice_energies,
    endpoint_singleton_slice_energies,
    middle_quintic_incidence_bound,
)
from repaired_open_profile_budget import coarse_open_completion_coefficients
from separated_balanced_endpoint_slice_contraction import (
    separated_balanced_coefficient,
    separated_balanced_orbit_entries,
)
from whole_cubic_middle_pair_contraction import (
    whole_cubic_middle_pair_coefficient,
    whole_cubic_middle_pair_orbit_entries,
)
from whole_cubic_quintic_triple_contraction import (
    whole_cubic_quintic_triple_coefficient,
    whole_cubic_quintic_triple_orbit_entries,
)


DIMENSION = 1024
THRESHOLD = 1 / 3
PROFILE = (1, 3, 5, 1)
SPLIT = (0, 2, 3, 0)
TARGET_GATE = 0.0529177166679669

Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class DoubleEndpointCubicQuinticRowContraction:
    dimension: int
    order: int
    cubic_pair_energy: float
    quintic_triple_energy: float
    record_one_cubic_incidence: int
    record_one_quintic_incidence: int
    record_one_middle_maximum: float
    record_one_row_bound: float
    record_three_middle_maximum: float
    record_three_row_bound: float
    row_energy_bound: float
    coefficient: float
    provisional_coefficient: float
    acceptance_gate: float
    optimal_beta: float
    optimized_total: float
    threshold_slack: float
    next_unresolved_entries: tuple[ProfileSplit, ...]
    next_unresolved_contribution: float
    next_admissible_coefficient: float


def double_endpoint_cubic_quintic_orbit_entries() -> tuple[ProfileSplit, ...]:
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


def scalar_row_sector_bounds(
    order: int,
) -> tuple[float, float, int, int, float, float, float, float]:
    """Return endpoint slices and the two scalar-row sector bounds."""

    q = order
    if q < 8 or q & (q - 1):
        raise ValueError(("power-of-two order at least eight required", q))
    dimension = q * q
    cubic_pair_energy = endpoint_singleton_slice_energies(q)[2]
    quintic_triple_energy = endpoint_quintic_singleton_slice_energies(q)[3]

    # In record one, the cubic is an L-shape.  At most 2(q-1) extensions
    # contain a fixed pair, each with endpoint amplitude 1/(q-1).  The
    # quintic is record one on both sides, with the exact fixed-triple
    # incidence below and endpoint amplitude at most one.
    cubic_incidence = 2 * (q - 1)
    quintic_incidence = middle_quintic_incidence_bound(q, 3)
    record_one_middle_maximum = (q + 2) / (q * (q - 1) * (q - 2))
    record_one_row = (
        cubic_incidence
        / (q - 1) ** 2
        * quintic_incidence
        * record_one_middle_maximum**2
    )

    # In record three, the unrestricted moment is at most 1/C(q,3) once
    # q>=8: the one-even-group correction is at most 3/(q-3)<=1.
    record_three_middle_maximum = 1 / comb(q, 3)
    record_three_row = (
        dimension**2
        * cubic_pair_energy
        * quintic_triple_energy
        * record_three_middle_maximum**2
    )
    return (
        cubic_pair_energy,
        quintic_triple_energy,
        cubic_incidence,
        quintic_incidence,
        record_one_middle_maximum,
        record_one_row,
        record_three_middle_maximum,
        record_three_row,
    )


def double_endpoint_cubic_quintic_coefficient(
    dimension: int = DIMENSION,
) -> float:
    """Return the proved arbitrary-diagonal coefficient for the orbit."""

    q = int(round(sqrt(dimension)))
    if q * q != dimension:
        raise ValueError(("square dimension required", dimension))
    bounds = scalar_row_sector_bounds(q)
    return sqrt(bounds[5] + bounds[7])


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


def double_endpoint_cubic_quintic_contraction(
    dimension: int = DIMENSION,
) -> DoubleEndpointCubicQuinticRowContraction:
    """Insert the tenth theorem, rerank, and compute the next gate."""

    if dimension != DIMENSION:
        raise ValueError(("ledger calibrated at N=1024", dimension))
    q = int(round(sqrt(dimension)))
    (
        cubic_pair_energy,
        quintic_triple_energy,
        cubic_incidence,
        quintic_incidence,
        record_one_maximum,
        record_one_row,
        record_three_maximum,
        record_three_row,
    ) = scalar_row_sector_bounds(q)
    row_energy = record_one_row + record_three_row
    coefficient = sqrt(row_energy)

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
        (
            whole_cubic_middle_pair_orbit_entries(),
            whole_cubic_middle_pair_coefficient(),
        ),
        (double_endpoint_cubic_quintic_orbit_entries(), coefficient),
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
    return DoubleEndpointCubicQuinticRowContraction(
        dimension=dimension,
        order=q,
        cubic_pair_energy=cubic_pair_energy,
        quintic_triple_energy=quintic_triple_energy,
        record_one_cubic_incidence=cubic_incidence,
        record_one_quintic_incidence=quintic_incidence,
        record_one_middle_maximum=record_one_maximum,
        record_one_row_bound=record_one_row,
        record_three_middle_maximum=record_three_maximum,
        record_three_row_bound=record_three_row,
        row_energy_bound=row_energy,
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
    result = double_endpoint_cubic_quintic_contraction()
    print(
        "double-endpoint cubic-quintic row contraction: "
        f"N={result.dimension},"
        f"q={result.order},"
        f"cubic_pair={result.cubic_pair_energy:.15g},"
        f"quintic_triple={result.quintic_triple_energy:.15g},"
        f"record_one_row={result.record_one_row_bound:.15g},"
        f"record_three_row={result.record_three_row_bound:.15g},"
        f"row_energy={result.row_energy_bound:.15g},"
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
