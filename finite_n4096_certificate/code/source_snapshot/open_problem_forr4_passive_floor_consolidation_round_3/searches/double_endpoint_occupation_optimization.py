#!/usr/bin/env python3
"""Optimize the complete dose-six (3,1,1,3) fixed-split ledger.

The objective retains one shared distribution over the 210 four-block
occupation states.  It combines all 64 ket/bra occurrence cuts with the
best currently proved coefficient for that cut.  A tangent hyperplane to
the concave geometric-mean objective supplies an a posteriori upper
certificate for the numerical optimizer.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import comb, sqrt

import cvxpy as cp
import numpy as np

from mixed_endpoint_weighted_bound import refined_q32_certificate
from same_middle_weighted_bound import bound as same_middle_bound


DOSE = 6
ORDER = 32
PROFILE = (3, 1, 1, 3)
SAME_ALTERNATING_UPPER = 0.010905


@dataclass(frozen=True)
class OccupationCertificate:
    objective: float
    supporting_upper: float
    attenuated_upper: float
    support: tuple[tuple[float, tuple[int, ...]], ...]
    maximum_gradient_state: tuple[int, ...]
    contribution_groups: tuple[tuple[str, float], ...]
    leading_cuts: tuple[tuple[float, tuple[int, ...]], ...]


def occupation_states() -> list[tuple[int, ...]]:
    return [
        state
        for state in product(range(DOSE + 1), repeat=4)
        if sum(state) <= DOSE
    ]


def splits() -> list[tuple[int, ...]]:
    return list(product(*(range(degree + 1) for degree in PROFILE)))


def feature(state: tuple[int, ...], split: tuple[int, ...]) -> int:
    result = 1
    for occupation, selected in zip(state, split, strict=True):
        result *= comb(occupation, selected)
    return result


def endpoint_energies(order: int) -> tuple[float, ...]:
    return (
        (order * order + 2) / 6,
        (order * order + 2) / (2 * order * order),
        (order * order - 2 * order + 2)
        / (order * order * (order - 1)),
        1 / (order * order),
    )


def old_endpoint_coefficient(left: int, right: int) -> float:
    energies = endpoint_energies(ORDER)
    value = min(
        sqrt(energies[left] * energies[right]),
        sqrt(energies[3 - left] * energies[3 - right]),
    )
    if {left, right} == {0, 3}:
        value = min(value, 1 / ORDER)
    return value


def block_coherent_coefficient(split: tuple[int, ...]) -> float:
    """Complete three-link cut table when every block is unsplit."""

    mask = frozenset(
        block
        for block, (degree, selected) in enumerate(
            zip(PROFILE, split, strict=True)
        )
        if selected == degree
    )
    if any(
        selected not in (0, degree)
        for degree, selected in zip(PROFILE, split, strict=True)
    ):
        raise ValueError(("not a block-coherent split", split))
    complement = frozenset(range(4)) - mask
    canonical = min(mask, complement, key=lambda value: (len(value), sorted(value)))
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


def coefficient(split: tuple[int, ...]) -> float:
    left, first_middle, second_middle, right = split
    balanced_endpoints = left in (1, 2) and right in (1, 2)
    same_endpoint_orientation = left == right
    same_middle_side = first_middle == second_middle

    if left in (0, 3) and right in (0, 3):
        return block_coherent_coefficient(split)

    if same_middle_side:
        if not balanced_endpoints:
            # The both-unsplit case was handled by the whole-block table.
            # The only remaining possibility has one balanced endpoint and
            # one unsplit endpoint.
            return float(same_middle_bound(ORDER).hybrid_upper)
        result = same_middle_bound(ORDER)
        if same_endpoint_orientation:
            return float(result.equal_endpoint_orientation)
        return float(result.mixed_endpoint_orientation)

    if balanced_endpoints:
        if same_endpoint_orientation:
            return SAME_ALTERNATING_UPPER
        return float(refined_q32_certificate().advertised_upper)
    return old_endpoint_coefficient(left, right)


def cut_group(split: tuple[int, ...]) -> str:
    left, first_middle, second_middle, right = split
    endpoint = (left, right)
    complement = (3 - left, 3 - right)
    canonical_endpoint = min(endpoint, complement)
    middle = "same_middle" if first_middle == second_middle else "alternating"
    if left in (1, 2) and right in (1, 2):
        orientation = "equal" if left == right else "mixed"
        return f"balanced_{orientation}_{middle}"
    return f"endpoint_{canonical_endpoint[0]}{canonical_endpoint[1]}_{middle}"


def certificate() -> OccupationCertificate:
    states = occupation_states()
    cut_list = splits()
    state_features = np.array(
        [
            [feature(state, split) for state in states]
            for split in cut_list
        ],
        dtype=float,
    )
    cut_coefficients = np.array(
        [coefficient(split) for split in cut_list], dtype=float
    )
    complement_indices = np.array(
        [
            cut_list.index(
                tuple(
                    degree - selected
                    for degree, selected in zip(
                        PROFILE, split, strict=True
                    )
                )
            )
            for split in cut_list
        ],
        dtype=int,
    )
    active_indices = [
        index
        for index, complement in enumerate(complement_indices)
        if np.any(state_features[index])
        and np.any(state_features[complement])
    ]

    weights = cp.Variable(len(states), nonneg=True)
    moments = state_features @ weights
    terms = [
        cut_coefficients[index]
        * cp.geo_mean(
            cp.hstack(
                [moments[index], moments[complement_indices[index]]]
            )
        )
        for index in active_indices
    ]
    problem = cp.Problem(cp.Maximize(cp.sum(terms)), [cp.sum(weights) == 1])
    objective = problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-9,
        tol_gap_rel=1e-9,
        tol_feas=1e-9,
        max_iter=1000,
    )
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError(problem.status)

    candidate = np.maximum(np.asarray(weights.value, dtype=float), 0)
    candidate /= candidate.sum()
    # Move an immaterial distance into the simplex interior so every
    # feasible cut moment has a finite supporting derivative.
    interior_mass = 1e-9
    candidate = (1 - interior_mass) * candidate + interior_mass / len(states)
    candidate_moments = state_features @ candidate
    if any(
        candidate_moments[index] <= 0
        or candidate_moments[complement_indices[index]] <= 0
        for index in active_indices
    ):
        raise AssertionError("candidate must have all cut moments positive")
    gradient = np.zeros(len(states))
    for index in active_indices:
        complement = complement_indices[index]
        ratio = sqrt(
            candidate_moments[complement]
            / candidate_moments[index]
        )
        gradient += cut_coefficients[index] / 2 * (
            ratio * state_features[index]
            + state_features[complement] / ratio
        )
    # The tangent inequality is exact; this allowance is deliberately much
    # larger than the observed floating-point variation in the moments and
    # gradients.  A formal interval certificate remains a separate cleanup.
    supporting_upper = float(np.max(gradient) + 1e-6)
    direct_objective = float(
        np.sum(
            cut_coefficients[active_indices]
            * np.sqrt(
                candidate_moments[active_indices]
                * candidate_moments[
                    complement_indices[active_indices]
                ]
            )
        )
    )
    if direct_objective > supporting_upper:
        raise AssertionError((direct_objective, supporting_upper))
    if abs(float(objective) - direct_objective) > 2e-7:
        raise AssertionError(("solver objective", objective, direct_objective))

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
    maximum_state = states[int(np.argmax(gradient))]
    cut_contributions = {
        cut_list[index]: float(
            cut_coefficients[index]
            * sqrt(
                candidate_moments[index]
                * candidate_moments[complement_indices[index]]
            )
        )
        for index in active_indices
    }
    grouped: dict[str, float] = {}
    for split, value in cut_contributions.items():
        group = cut_group(split)
        grouped[group] = grouped.get(group, 0.0) + value
    return OccupationCertificate(
        objective=direct_objective,
        supporting_upper=supporting_upper,
        attenuated_upper=(5 / 6) ** 8 * supporting_upper,
        support=support,
        maximum_gradient_state=maximum_state,
        contribution_groups=tuple(
            sorted(grouped.items(), key=lambda item: item[1], reverse=True)
        ),
        leading_cuts=tuple(
            sorted(
                ((value, split) for split, value in cut_contributions.items()),
                reverse=True,
            )[:12]
        ),
    )


def main() -> None:
    result = certificate()
    print(
        "complete double-endpoint occupation optimization: "
        f"objective={result.objective:.15g},"
        f"supporting_upper={result.supporting_upper:.15g},"
        f"attenuated_upper={result.attenuated_upper:.15g},"
        f"maximum_gradient_state={result.maximum_gradient_state}"
    )
    print(
        "support="
        + ";".join(
            f"{state}:{weight:.12g}" for weight, state in result.support
        )
    )
    print(
        "contribution_groups="
        + ";".join(
            f"{group}:{value:.12g}"
            for group, value in result.contribution_groups
        )
    )
    print(
        "leading_cuts="
        + ";".join(
            f"{split}:{value:.12g}"
            for value, split in result.leading_cuts
        )
    )


if __name__ == "__main__":
    main()
