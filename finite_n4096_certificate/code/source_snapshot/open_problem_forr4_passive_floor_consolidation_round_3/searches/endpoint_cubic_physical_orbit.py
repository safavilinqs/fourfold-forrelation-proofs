#!/usr/bin/env python3
"""Direct physical-orbit diagnostic for the two endpoint-cubic profiles.

The fixed-mask ledger optimizes each Schur flattening after forgetting that
all masks arise from one passive diagonal law.  This script keeps that law:
it first optimizes the current q=32 occupation relaxation for the two
endpoint-cubic profiles, then realizes the resulting occupation mixture as
the coordinate-uniform passive law at q=2.  On that actual support basis it
builds the complete Hermitian endpoint-cubic Schur kernel, with both
(3,1,1,1) and (1,1,1,3) entries present simultaneously, and computes its
trace norm.

This is a falsification/structure diagnostic, not a q=32 certificate.  Its
purpose is to measure the loss from fixed-mask triangle inequality and to
identify a tractable compound theorem target.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
from math import comb, prod, sqrt

import cvxpy as cp
import numpy as np

from double_endpoint_occupation_optimization import occupation_states
from single_cubic_weighted_bound import bound as endpoint_cubic_bound


DOSE = 6
TARGET_ORDER = 32
DIRECT_ORDER = 2
ENDPOINT_PROFILES = ((3, 1, 1, 1), (1, 1, 1, 3))


@dataclass(frozen=True)
class OccupationCandidate:
    objective: float
    weights: tuple[tuple[tuple[int, ...], float], ...]


@dataclass(frozen=True)
class PhysicalOrbitResult:
    dimension: int
    nonzero_entries: int
    separate_nuclear: tuple[float, float]
    joint_nuclear: float
    q2_fixed_mask_upper: float
    occupation: OccupationCandidate


def state_feature(state: tuple[int, ...], split: tuple[int, ...]) -> int:
    return prod(
        comb(occupation, selected)
        for occupation, selected in zip(state, split, strict=True)
    )


def splits(profile: tuple[int, ...]) -> list[tuple[int, ...]]:
    return list(product(*(range(degree + 1) for degree in profile)))


def path_cut_coefficient(order: int, mask: frozenset[int]) -> float:
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
    return order ** (-exponent)


def coefficient(
    order: int, profile: tuple[int, ...], split: tuple[int, ...]
) -> float:
    decorated_block = 0 if profile[0] == 3 else 3
    if split[decorated_block] in (0, 3):
        if all(
            selected in (0, degree)
            for degree, selected in zip(profile, split, strict=True)
        ):
            mask = frozenset(
                block
                for block, (degree, selected) in enumerate(
                    zip(profile, split, strict=True)
                )
                if selected == degree
            )
            return path_cut_coefficient(order, mask)
        raise AssertionError((profile, split))
    singleton_count = sum(
        selected
        for block, selected in enumerate(split)
        if block != decorated_block
    )
    result = endpoint_cubic_bound(order)
    if singleton_count in (0, 3):
        return float(result.extreme_singletons)
    return float(result.balanced_singletons)


def optimize_occupation(order: int) -> OccupationCandidate:
    states = occupation_states()
    rows: list[np.ndarray] = []
    complement_rows: list[np.ndarray] = []
    constants: list[float] = []
    for profile in ENDPOINT_PROFILES:
        split_list = splits(profile)
        index = {split: position for position, split in enumerate(split_list)}
        features = np.asarray(
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
            if not np.any(features[position]) or not np.any(
                features[index[complement]]
            ):
                continue
            rows.append(features[position])
            complement_rows.append(features[index[complement]])
            constants.append(coefficient(order, profile, split))

    left = np.asarray(rows)
    right = np.asarray(complement_rows)
    constants_array = np.asarray(constants)
    weights = cp.Variable(len(states), nonneg=True)
    left_moments = left @ weights
    right_moments = right @ weights
    objective = cp.sum(
        [
            constants_array[index]
            * cp.geo_mean(
                cp.hstack([left_moments[index], right_moments[index]])
            )
            for index in range(len(constants_array))
        ]
    )
    problem = cp.Problem(cp.Maximize(objective), [cp.sum(weights) == 1])
    value = problem.solve(
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
    support = tuple(
        (state, float(weight))
        for state, weight in zip(states, candidate, strict=True)
        if weight > 2e-7
    )
    retained_mass = sum(weight for _, weight in support)
    support = tuple((state, weight / retained_mass) for state, weight in support)
    return OccupationCandidate(float(value), support)


def subsets_by_size(dimension: int) -> list[list[int]]:
    result: list[list[int]] = []
    for size in range(DOSE + 1):
        if size > dimension:
            result.append([])
            continue
        result.append(
            [
                sum(1 << coordinate for coordinate in chosen)
                for chosen in combinations(range(dimension), size)
            ]
        )
    return result


def sylvester(order: int) -> np.ndarray:
    result = np.asarray([[1.0]])
    while len(result) < order:
        result = np.block([[result, result], [result, -result]])
    return result / sqrt(order)


def endpoint_link(
    order: int, endpoint_support: int, singleton: int
) -> float:
    dimension = order * order
    coordinates = [
        coordinate
        for coordinate in range(dimension)
        if (endpoint_support >> coordinate) & 1
    ]
    if len(coordinates) != 3 or singleton.bit_count() != 1:
        raise ValueError((endpoint_support, singleton))
    hidden = [coordinate % order for coordinate in coordinates]
    distinct = len(set(hidden))
    if distinct == 1:
        weight = 1.0
    elif distinct == 2:
        weight = -1 / (order - 1)
    else:
        return 0.0
    xor_coordinate = coordinates[0] ^ coordinates[1] ^ coordinates[2]
    singleton_coordinate = singleton.bit_length() - 1
    return weight * sylvester(dimension)[xor_coordinate, singleton_coordinate]


def kernel_value(
    order: int, differences: tuple[int, ...], profile: tuple[int, ...]
) -> float:
    singletons = [value.bit_length() - 1 for value in differences]
    hadamard = sylvester(order * order)
    if profile[0] == 3:
        return (
            endpoint_link(order, differences[0], differences[1])
            * hadamard[singletons[1], singletons[2]]
            * hadamard[singletons[2], singletons[3]]
        )
    return (
        hadamard[singletons[0], singletons[1]]
        * hadamard[singletons[1], singletons[2]]
        * endpoint_link(order, differences[3], differences[2])
    )


def occupation_upper(
    order: int,
    law: tuple[tuple[tuple[int, ...], float], ...],
) -> float:
    total = 0.0
    for profile in ENDPOINT_PROFILES:
        for split in splits(profile):
            complement = tuple(
                degree - selected
                for degree, selected in zip(profile, split, strict=True)
            )
            left = sum(
                weight * state_feature(state, split) for state, weight in law
            )
            right = sum(
                weight * state_feature(state, complement)
                for state, weight in law
            )
            if left and right:
                total += coefficient(order, profile, split) * sqrt(left * right)
    return total


def direct_physical_orbit() -> PhysicalOrbitResult:
    occupation = optimize_occupation(TARGET_ORDER)
    dimension = DIRECT_ORDER * DIRECT_ORDER
    by_size = subsets_by_size(dimension)
    basis: list[tuple[int, ...]] = []
    probabilities: list[float] = []
    for state, weight in occupation.weights:
        if any(size > dimension for size in state):
            continue
        orbit = list(product(*(by_size[size] for size in state)))
        if not orbit:
            continue
        basis.extend(orbit)
        probabilities.extend([weight / len(orbit)] * len(orbit))
    total_probability = sum(probabilities)
    probabilities = [value / total_probability for value in probabilities]
    law = tuple(
        (state, weight / total_probability)
        for state, weight in occupation.weights
        if all(size <= dimension for size in state)
    )

    matrices = [np.zeros((len(basis), len(basis))) for _ in ENDPOINT_PROFILES]
    nonzero = 0
    for row in range(len(basis)):
        for column in range(row + 1, len(basis)):
            differences = tuple(
                left ^ right
                for left, right in zip(basis[row], basis[column], strict=True)
            )
            sizes = tuple(value.bit_count() for value in differences)
            if sizes not in ENDPOINT_PROFILES:
                continue
            profile_index = ENDPOINT_PROFILES.index(sizes)
            value = kernel_value(DIRECT_ORDER, differences, sizes)
            if not value:
                continue
            weighted = sqrt(probabilities[row] * probabilities[column]) * value
            matrices[profile_index][row, column] = weighted
            matrices[profile_index][column, row] = weighted
            nonzero += 2

    separate = tuple(
        float(np.abs(np.linalg.eigvalsh(matrix)).sum()) for matrix in matrices
    )
    joint = float(np.abs(np.linalg.eigvalsh(matrices[0] + matrices[1])).sum())
    return PhysicalOrbitResult(
        dimension=len(basis),
        nonzero_entries=nonzero,
        separate_nuclear=(separate[0], separate[1]),
        joint_nuclear=joint,
        q2_fixed_mask_upper=occupation_upper(DIRECT_ORDER, law),
        occupation=occupation,
    )


def main() -> None:
    result = direct_physical_orbit()
    print(
        "endpoint-cubic physical-orbit diagnostic: "
        f"q32_relaxation={result.occupation.objective:.12g},"
        f"q2_dimension={result.dimension},"
        f"q2_nonzero={result.nonzero_entries},"
        f"q2_fixed_mask_upper={result.q2_fixed_mask_upper:.12g},"
        f"q2_separate={result.separate_nuclear[0]:.12g},"
        f"{result.separate_nuclear[1]:.12g},"
        f"q2_joint={result.joint_nuclear:.12g},"
        f"joint_to_fixed={result.joint_nuclear/result.q2_fixed_mask_upper:.12g}"
    )
    print(
        "occupation="
        + ";".join(
            f"{state}:{weight:.12g}"
            for state, weight in result.occupation.weights
        )
    )


if __name__ == "__main__":
    main()
