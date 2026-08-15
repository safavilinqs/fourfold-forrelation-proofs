#!/usr/bin/env python3
"""Shared-law contraction for the fourth balanced high-sector orbit.

The target is profile ``(3,1,1,5)`` and split ``(1,1,1,2)``.  The two
singleton blocks lie on the same side of the cut.  A symmetry-twirled exact
coefficient for the split cubic endpoint is dressed by a completed-Gram
minus overlap factor for the split quintic endpoint.  The internal singleton
Hadamard link supplies ``1/q``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import sqrt

import numpy as np
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
    leading_balanced_coefficient,
    leading_balanced_orbit_entries,
)
from occupation_compatible_sector_optimization import (
    certificate,
    endpoint_quintic_singleton_slice_energies,
    endpoint_singleton_slice_energies,
)
from opposite_endpoint_orbit_scan import (
    cubic_response,
    cubic_response_summary,
    quintic_weight,
    support_xor,
    walsh_transform,
)
from repaired_open_profile_budget import coarse_open_completion_coefficients
from separated_balanced_endpoint_slice_contraction import (
    separated_balanced_coefficient,
    separated_balanced_orbit_entries,
)


DIMENSION = 1024
THRESHOLD = 1 / 3
PROFILE = (3, 1, 1, 5)
SPLIT = (1, 1, 1, 2)

Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class CubicTwirl:
    order: int
    vertical_spectrum_numerator: int
    horizontal_spectrum_numerator: int
    general_spectrum_numerator: int
    vertical_contribution: float
    horizontal_contribution: float
    general_contribution: float
    coefficient: float


@dataclass(frozen=True)
class InternalSingletonContraction:
    dimension: int
    order: int
    cubic_twirled_coefficient: float
    quintic_overlap_factor: float
    coefficient: float
    provisional_coefficient: float
    acceptance_gate: float
    simple_slice_coefficient: float
    simple_slice_optimized_total: float
    simple_slice_threshold_overshoot: float
    optimal_beta: float
    optimized_total: float
    threshold_slack: float
    next_unresolved_entries: tuple[ProfileSplit, ...]
    next_unresolved_contribution: float
    next_admissible_coefficient: float
    vertical_mixture_diagnostic: float


def internal_singleton_orbit_entries() -> tuple[ProfileSplit, ...]:
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


def _integer_cubic_spectrum_numerator(order: int, difference: int) -> int:
    """Return the exact allowed Walsh l1 numerator for one pair orbit."""

    dimension = order * order
    response = cubic_response(order, difference, False)
    integer_response = np.rint(response * (order - 1)).astype(np.int64)
    spectrum = np.rint(walsh_transform(integer_response)).astype(np.int64)
    allowed = np.asarray(
        [
            int(difference & frequency).bit_count() % 2 == 0
            for frequency in range(dimension)
        ]
    )
    return int(np.abs(spectrum[allowed]).sum())


def cubic_twirled_coefficient(order: int = 32) -> CubicTwirl:
    """Exact arbitrary-law coefficient for the split cubic endpoint.

    Root fidelity is jointly concave in the row and column diagonal laws.
    Twirling therefore makes the row law uniform and leaves three pair
    orbits: vertical, horizontal, and general.  Their orthogonal Walsh blocks
    reduce the remaining optimization to Cauchy--Schwarz in three weights.
    """

    if order < 4 or order & (order - 1):
        raise ValueError(("power-of-two order at least four required", order))
    q = order
    dimension = q * q
    representatives = (q, 1, q + 1)
    difference_counts = (q - 1, q - 1, (q - 1) ** 2)
    numerators = tuple(
        _integer_cubic_spectrum_numerator(q, difference)
        for difference in representatives
    )
    contributions = tuple(
        numerator * sqrt(count)
        / ((q - 1) * dimension * sqrt(dimension))
        for numerator, count in zip(
            numerators, difference_counts, strict=True
        )
    )
    return CubicTwirl(
        order=q,
        vertical_spectrum_numerator=numerators[0],
        horizontal_spectrum_numerator=numerators[1],
        general_spectrum_numerator=numerators[2],
        vertical_contribution=contributions[0],
        horizontal_contribution=contributions[1],
        general_contribution=contributions[2],
        coefficient=sqrt(sum(value * value for value in contributions)),
    )


def internal_singleton_coefficient(dimension: int = DIMENSION) -> float:
    """Return the proved arbitrary-diagonal coefficient for the orbit."""

    order = int(round(sqrt(dimension)))
    if order * order != dimension:
        raise ValueError(("square dimension required", dimension))
    cubic = cubic_twirled_coefficient(order).coefficient
    return cubic * (1 + sqrt(2)) / order


def vertical_mixture_diagnostic(order: int = 32) -> float:
    """Exact uniform vertical-family coefficient; not an upper theorem."""

    if order < 4 or order & (order - 1):
        raise ValueError(("power-of-two order at least four required", order))
    q = order
    dimension = q * q
    vertical_pairs = tuple(
        tuple(row * q + column for row in rows)
        for column in range(q)
        for rows in combinations(range(q), 2)
    )
    cubic_nuclear = (
        (q - 1) * cubic_response_summary(q, q, False).nuclear
    )
    triples = tuple(
        tuple(row * q + column for row in rows)
        for column in range(q)
        for rows in combinations(range(q), 3)
    )
    zero_xor_triples = tuple(
        triple for triple in triples if support_xor(triple) == 0
    )
    quintic_block = np.zeros(
        (len(vertical_pairs), len(zero_xor_triples))
    )
    for pair_index, pair in enumerate(vertical_pairs):
        pair_set = set(pair)
        for triple_index, triple in enumerate(zero_xor_triples):
            if pair_set.intersection(triple):
                continue
            quintic_block[pair_index, triple_index] = quintic_weight(
                tuple(sorted(pair + triple)), q
            )
    quintic_nuclear = dimension * float(
        np.linalg.svd(quintic_block, compute_uv=False).sum()
    )
    cubic_normalized = cubic_nuclear / sqrt(
        dimension * dimension * len(vertical_pairs)
    )
    quintic_normalized = quintic_nuclear / sqrt(
        len(vertical_pairs) * dimension * len(triples)
    )
    return cubic_normalized * quintic_normalized / q


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


def internal_singleton_contraction(
    dimension: int = DIMENSION,
) -> InternalSingletonContraction:
    """Insert the fourth theorem, rerank, and compute the next gate."""

    if dimension != DIMENSION:
        raise ValueError(("ledger calibrated at N=1024", dimension))
    q = int(round(sqrt(dimension)))
    cubic = cubic_twirled_coefficient(q)
    coefficient = internal_singleton_coefficient(dimension)
    base_coefficients = coarse_open_completion_coefficients()
    coefficients = dict(base_coefficients)
    provisional = base_coefficients[(PROFILE, SPLIT)]
    for orbit, value in (
        (leading_balanced_orbit_entries(), leading_balanced_coefficient()),
        (adjacent_balanced_orbit_entries(), adjacent_balanced_coefficient()),
        (separated_balanced_orbit_entries(), separated_balanced_coefficient()),
        (internal_singleton_orbit_entries(), coefficient),
    ):
        for entry in orbit:
            coefficients[entry] = value

    cubic_slices = endpoint_singleton_slice_energies(q)
    quintic_slices = endpoint_quintic_singleton_slice_energies(q)
    simple_slice = sqrt(
        cubic_slices[1] * quintic_slices[2] / dimension
    )
    simple_coefficients = dict(coefficients)
    for entry in internal_singleton_orbit_entries():
        simple_coefficients[entry] = simple_slice
    _, simple_total = _optimize_total(simple_coefficients, dimension)

    optimal_beta, optimized_total = _optimize_total(coefficients, dimension)
    ledger = certificate(
        beta=optimal_beta,
        profile_split_coefficients=coefficients,
    )
    weights = dict(ledger.occupation_weights)
    proved_orbits = tuple(
        frozenset(orbit)
        for orbit in (
            leading_balanced_orbit_entries(),
            adjacent_balanced_orbit_entries(),
            separated_balanced_orbit_entries(),
            internal_singleton_orbit_entries(),
        )
    )
    unresolved = []
    for orbit in symmetry_orbits(dose_six_relevant_entries()):
        if frozenset(orbit) in proved_orbits:
            continue
        if not all(
            abs(base_coefficients[entry] - 1 / q) < 1e-14
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
        1 / q,
        0.2,
        xtol=1e-14,
    )
    return InternalSingletonContraction(
        dimension=dimension,
        order=q,
        cubic_twirled_coefficient=cubic.coefficient,
        quintic_overlap_factor=1 + sqrt(2),
        coefficient=coefficient,
        provisional_coefficient=provisional,
        acceptance_gate=0.0450405467777823,
        simple_slice_coefficient=simple_slice,
        simple_slice_optimized_total=simple_total,
        simple_slice_threshold_overshoot=simple_total - THRESHOLD,
        optimal_beta=optimal_beta,
        optimized_total=optimized_total,
        threshold_slack=THRESHOLD - optimized_total,
        next_unresolved_entries=next_entries,
        next_unresolved_contribution=next_contribution,
        next_admissible_coefficient=next_admissible,
        vertical_mixture_diagnostic=vertical_mixture_diagnostic(q),
    )


def main() -> None:
    cubic = cubic_twirled_coefficient()
    result = internal_singleton_contraction()
    print(
        "internal singleton shared-law contraction: "
        f"N={result.dimension},q={result.order},"
        f"cubic_numerators="
        f"{cubic.vertical_spectrum_numerator}/"
        f"{cubic.horizontal_spectrum_numerator}/"
        f"{cubic.general_spectrum_numerator},"
        f"cubic_twirled={result.cubic_twirled_coefficient:.15g},"
        f"quintic_overlap={result.quintic_overlap_factor:.15g},"
        f"coefficient={result.coefficient:.15g},"
        f"gate={result.acceptance_gate:.15g},"
        f"simple_slice={result.simple_slice_coefficient:.15g},"
        f"simple_slice_total={result.simple_slice_optimized_total:.15g},"
        f"simple_slice_overshoot="
        f"{result.simple_slice_threshold_overshoot:.15g},"
        f"optimal_beta={result.optimal_beta:.15g},"
        f"optimized_total={result.optimized_total:.15g},"
        f"threshold_slack={result.threshold_slack:.15g},"
        f"next_unresolved={result.next_unresolved_entries[0]},"
        f"next_contribution={result.next_unresolved_contribution:.15g},"
        f"next_admissible={result.next_admissible_coefficient:.15g},"
        f"vertical_mixture={result.vertical_mixture_diagnostic:.15g}"
    )


if __name__ == "__main__":
    main()
