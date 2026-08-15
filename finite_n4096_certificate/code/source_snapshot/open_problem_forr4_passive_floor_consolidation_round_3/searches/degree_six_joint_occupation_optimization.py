#!/usr/bin/env python3
"""Joint occupation diagnostic for degree six plus the double endpoint.

This is the first A3 ledger.  It retains one law on all 210 dose-six
occupation states, inserts the refined (3,1,1,3) cut table, and adds the four
single-cubic degree-six profiles with their currently proved fixed-cut
constants.  A tangent hyperplane certifies the displayed numerical upper.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import comb

import cvxpy as cp
import numpy as np

from double_endpoint_occupation_optimization import (
    DOSE,
    PROFILE as DOUBLE_ENDPOINT_PROFILE,
    coefficient as double_endpoint_coefficient,
    occupation_states,
)
from single_cubic_weighted_bound import bound as endpoint_cubic_bound


ORDER = 32
BETA = 5 / 6
DEGREE_SIX_PROFILES = (
    (3, 1, 1, 1),
    (1, 3, 1, 1),
    (1, 1, 3, 1),
    (1, 1, 1, 3),
)


@dataclass(frozen=True)
class JointDegreeSixCertificate:
    objective: float
    supporting_upper: float
    double_endpoint_at_candidate: float
    degree_six_at_candidate: float
    support: tuple[tuple[float, tuple[int, ...]], ...]
    maximum_gradient_state: tuple[int, ...]
    profile_contributions: tuple[tuple[tuple[int, ...], float], ...]
    leading_terms: tuple[
        tuple[float, tuple[int, ...], tuple[int, ...]], ...
    ]


def path_cut_coefficient(mask: frozenset[int]) -> float:
    complement = frozenset(range(4)) - mask
    canonical = min(
        mask, complement, key=lambda value: (len(value), sorted(value))
    )
    if not canonical:
        exponent = 3
    elif len(canonical) == 1:
        exponent = 2
    else:
        exponent = 2 if canonical in (
            frozenset({0, 1}),
            frozenset({2, 3}),
        ) else 1
    return ORDER ** (-exponent)


def middle_cubic_path_coefficient(
    decorated_block: int, mask: frozenset[int]
) -> float:
    """Whole-block table with two L-shape links and one Hadamard link."""

    if decorated_block not in (1, 2):
        raise ValueError(decorated_block)
    if decorated_block == 2:
        mask = frozenset(3 - block for block in mask)
    q = ORDER
    n = q - 1
    base = {
        frozenset(): 1 / (q**3 * n**2),
        frozenset({0}): 1 / (q**2 * n),
        frozenset({1}): 1 / (q**2 * n),
        frozenset({2}): 1 / (q**2 * n**2),
        frozenset({3}): 1 / (q**2 * n**2),
        frozenset({0, 1}): 1 / (q**2 * n),
        frozenset({1, 2}): 1 / (q * n),
        frozenset({0, 2}): 1 / (q * n),
    }
    if mask in base:
        return base[mask]
    complement = frozenset(range(4)) - mask
    if complement not in base:
        raise AssertionError(("middle cubic whole-block mask", mask))
    return base[complement]


def degree_six_coefficient(
    profile: tuple[int, ...], split: tuple[int, ...]
) -> float:
    block_coherent = all(
        selected in (0, degree)
        for degree, selected in zip(profile, split, strict=True)
    )
    if block_coherent:
        mask = frozenset(
            block
            for block, (degree, selected) in enumerate(
                zip(profile, split, strict=True)
            )
            if selected == degree
        )
        decorated_block = next(
            block for block, degree in enumerate(profile) if degree == 3
        )
        if decorated_block in (1, 2):
            return middle_cubic_path_coefficient(decorated_block, mask)
        return path_cut_coefficient(mask)
    decorated_block = next(
        block for block, degree in enumerate(profile) if degree == 3
    )
    if decorated_block in (0, 3):
        singleton_count = sum(
            selected
            for block, selected in enumerate(split)
            if block != decorated_block
        )
        result = endpoint_cubic_bound(ORDER)
        if singleton_count in (0, 3):
            return float(result.extreme_singletons)
        return float(result.balanced_singletons)
    # For a cubic middle block, the two singleton records force an L-shape.
    # Every entry is 1/[q^3(q-1)^2], and rank--Frobenius over six marks gives
    # the arbitrary-diagonal coefficient 1/(q-1)^2.
    return 1 / (ORDER - 1) ** 2


def profile_splits(profile: tuple[int, ...]) -> list[tuple[int, ...]]:
    return list(product(*(range(degree + 1) for degree in profile)))


def state_feature(
    state: tuple[int, ...], split: tuple[int, ...]
) -> int:
    result = 1
    for occupation, selected in zip(state, split, strict=True):
        result *= comb(occupation, selected)
    return result


def certificate() -> JointDegreeSixCertificate:
    states = occupation_states()
    term_rows: list[np.ndarray] = []
    complement_rows: list[np.ndarray] = []
    term_coefficients: list[float] = []
    term_degrees: list[int] = []
    term_profiles: list[tuple[int, ...]] = []
    term_splits: list[tuple[int, ...]] = []

    profiles = (DOUBLE_ENDPOINT_PROFILE,) + DEGREE_SIX_PROFILES
    for profile in profiles:
        split_list = profile_splits(profile)
        index = {split: position for position, split in enumerate(split_list)}
        features = np.array(
            [
                [state_feature(state, split) for state in states]
                for split in split_list
            ],
            dtype=float,
        )
        for position, split in enumerate(split_list):
            complement = tuple(
                degree - selected
                for degree, selected in zip(profile, split, strict=True)
            )
            complement_position = index[complement]
            if not np.any(features[position]) or not np.any(
                features[complement_position]
            ):
                continue
            coefficient = (
                double_endpoint_coefficient(split)
                if profile == DOUBLE_ENDPOINT_PROFILE
                else degree_six_coefficient(profile, split)
            )
            term_rows.append(features[position])
            complement_rows.append(features[complement_position])
            term_coefficients.append(coefficient * BETA ** sum(profile))
            term_degrees.append(sum(profile))
            term_profiles.append(profile)
            term_splits.append(split)

    left_features = np.array(term_rows)
    right_features = np.array(complement_rows)
    coefficients = np.array(term_coefficients)
    degrees = np.array(term_degrees)
    weights = cp.Variable(len(states), nonneg=True)
    left_moments = left_features @ weights
    right_moments = right_features @ weights
    objective_terms = [
        coefficients[index]
        * cp.geo_mean(
            cp.hstack([left_moments[index], right_moments[index]])
        )
        for index in range(len(coefficients))
    ]
    problem = cp.Problem(
        cp.Maximize(cp.sum(objective_terms)), [cp.sum(weights) == 1]
    )
    solver_objective = problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-8,
        tol_gap_rel=1e-8,
        tol_feas=1e-8,
        max_iter=1000,
    )
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError(problem.status)

    candidate = np.maximum(np.asarray(weights.value, dtype=float), 0)
    candidate /= candidate.sum()
    interior_mass = 1e-9
    candidate = (1 - interior_mass) * candidate + interior_mass / len(states)
    left = left_features @ candidate
    right = right_features @ candidate
    if np.any(left <= 0) or np.any(right <= 0):
        raise AssertionError("interior candidate moments")
    values = coefficients * np.sqrt(left * right)
    direct_objective = float(values.sum())
    if abs(float(solver_objective) - direct_objective) > 3e-6:
        raise AssertionError((solver_objective, direct_objective))

    gradient = np.zeros(len(states))
    for index in range(len(coefficients)):
        ratio = np.sqrt(right[index] / left[index])
        gradient += coefficients[index] / 2 * (
            ratio * left_features[index]
            + right_features[index] / ratio
        )
    supporting_upper = float(np.max(gradient) + 1e-6)
    if direct_objective > supporting_upper:
        raise AssertionError((direct_objective, supporting_upper))

    support = tuple(
        sorted(
            (
                (float(weight), state)
                for weight, state in zip(candidate, states, strict=True)
                if weight > 2e-7
            ),
            reverse=True,
        )
    )
    profile_totals: dict[tuple[int, ...], float] = {}
    for profile, value in zip(term_profiles, values, strict=True):
        profile_totals[profile] = profile_totals.get(profile, 0.0) + float(
            value
        )
    return JointDegreeSixCertificate(
        objective=direct_objective,
        supporting_upper=supporting_upper,
        double_endpoint_at_candidate=float(values[degrees == 8].sum()),
        degree_six_at_candidate=float(values[degrees == 6].sum()),
        support=support,
        maximum_gradient_state=states[int(np.argmax(gradient))],
        profile_contributions=tuple(
            sorted(
                profile_totals.items(),
                key=lambda item: item[1],
                reverse=True,
            )
        ),
        leading_terms=tuple(
            sorted(
                (
                    (float(value), profile, split)
                    for value, profile, split in zip(
                        values, term_profiles, term_splits, strict=True
                    )
                ),
                reverse=True,
            )[:16]
        ),
    )


def main() -> None:
    result = certificate()
    print(
        "degree-six plus double-endpoint occupation diagnostic: "
        f"objective={result.objective:.15g},"
        f"supporting_upper={result.supporting_upper:.15g},"
        f"double_endpoint={result.double_endpoint_at_candidate:.15g},"
        f"degree_six={result.degree_six_at_candidate:.15g},"
        f"available_margin={0.160358131958:.15g},"
        f"overshoot={result.supporting_upper-0.160358131958:.15g},"
        f"maximum_gradient_state={result.maximum_gradient_state}"
    )
    print(
        "support="
        + ";".join(
            f"{state}:{weight:.12g}" for weight, state in result.support
        )
    )
    print(
        "profile_contributions="
        + ";".join(
            f"{profile}:{value:.12g}"
            for profile, value in result.profile_contributions
        )
    )
    print(
        "leading_terms="
        + ";".join(
            f"{profile}/{split}:{value:.12g}"
            for value, profile, split in result.leading_terms
        )
    )


if __name__ == "__main__":
    main()
