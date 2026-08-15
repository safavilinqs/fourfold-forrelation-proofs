#!/usr/bin/env python3
"""Independent physical audit of the twelve dual-endpoint entries."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))
sys.path.append(str(ROOT / "tests"))

from q64_dual_endpoint_schur_insertion import (  # noqa: E402
    cubic_fixed_pair_energy,
    dual_endpoint_entries,
    quintic_fixed_triple_energy,
)
from q64_recovered_cubic_quintic_independent_audit import (  # noqa: E402
    DirectQ4Plant,
    direct_permutation_moment,
    transpose_support,
)
from signed_permutation_link_moment import moment  # noqa: E402


SEED = 2026071712
Q4 = 4


def exact_q4_endpoint_energies(plant: DirectQ4Plant) -> tuple[Fraction, Fraction]:
    """Enumerate the two physical simple-support endpoint slices."""

    cubic = plant.moments(3, 1).astype(np.int64)
    quintic = plant.moments(1, 5).astype(np.int64)
    cubic_maximum = 0
    for pair in combinations(range(Q4 * Q4), 2):
        pair_set = set(pair)
        extensions = [
            plant.support_index[3][tuple(sorted(pair + (cell,)))]
            for cell in range(Q4 * Q4)
            if cell not in pair_set
        ]
        cubic_maximum = max(
            cubic_maximum,
            max(
                int(np.square(cubic[extensions, singleton]).sum())
                for singleton in range(Q4 * Q4)
            ),
        )

    quintic_maximum = 0
    for triple in combinations(range(Q4 * Q4), 3):
        available = tuple(cell for cell in range(Q4 * Q4) if cell not in triple)
        extensions = [
            plant.support_index[5][tuple(sorted(triple + pair))]
            for pair in combinations(available, 2)
        ]
        quintic_maximum = max(
            quintic_maximum,
            max(
                int(np.square(quintic[singleton, extensions]).sum())
                for singleton in range(Q4 * Q4)
            ),
        )
    scale = plant.group_size**2
    return Fraction(cubic_maximum, scale), Fraction(quintic_maximum, scale)


def partial_feature(values: np.ndarray, support: tuple[int, ...]) -> np.ndarray:
    if not support:
        return np.ones(values.shape[0], dtype=np.int16)
    return np.prod(values[:, support], axis=1, dtype=np.int16)


def numerator(
    plant: DirectQ4Plant, left: tuple[int, ...], right: tuple[int, ...]
) -> int:
    return int(
        plant.moments(len(left), len(right))[
            plant.support_index[len(left)][left],
            plant.support_index[len(right)][right],
        ]
    )


def lifted_gram_checks(plant: DirectQ4Plant) -> int:
    """Check the two nontrivial completed-link Gram lifts exactly."""

    rng = np.random.default_rng(SEED)
    checks = 0
    for _ in range(96):
        quintic = plant.supports[5][int(rng.integers(len(plant.supports[5])))]
        left_pair = tuple(sorted(rng.choice(quintic, size=2, replace=False)))
        right_triple = tuple(cell for cell in quintic if cell not in left_pair)
        cubic = plant.supports[3][int(rng.integers(len(plant.supports[3])))]
        lifted = int(
            (
                partial_feature(plant.left_values, left_pair)
                * partial_feature(plant.right_values, cubic)
            )
            @ partial_feature(plant.left_values, right_triple)
        )
        if lifted != numerator(plant, quintic, cubic):
            raise AssertionError(("quintic completed-link Gram lift", quintic, cubic))
        checks += 1

        split_cubic = plant.supports[3][
            int(rng.integers(len(plant.supports[3])))
        ]
        row_pair = tuple(sorted(rng.choice(split_cubic, size=2, replace=False)))
        column_cell = tuple(cell for cell in split_cubic if cell not in row_pair)
        whole_cubic = plant.supports[3][
            int(rng.integers(len(plant.supports[3])))
        ]
        lifted = int(
            partial_feature(plant.right_values, row_pair)
            @ (
                partial_feature(plant.left_values, whole_cubic)
                * partial_feature(plant.right_values, column_cell)
            )
        )
        if lifted != numerator(plant, whole_cubic, split_cubic):
            raise AssertionError(
                ("cubic completed-link Gram lift", whole_cubic, split_cubic)
            )
        checks += 1
    return checks


def random_partial(
    rng: np.random.Generator, size: int, dimension: int = Q4 * Q4
) -> tuple[int, ...]:
    return tuple(sorted(int(value) for value in rng.choice(dimension, size=size, replace=False)))


def occurrence_submatrix(
    plant: DirectQ4Plant,
    entry: tuple[tuple[int, ...], tuple[int, ...]],
    rng: np.random.Generator,
    row_count: int = 28,
    column_count: int = 36,
) -> np.ndarray:
    profile, split = entry
    rows = tuple(
        tuple(random_partial(rng, selected) for selected in split)
        for _ in range(row_count)
    )
    columns = tuple(
        tuple(
            random_partial(rng, degree - selected)
            for degree, selected in zip(profile, split, strict=True)
        )
        for _ in range(column_count)
    )
    matrix = np.zeros((row_count, column_count))
    for row_index, row in enumerate(rows):
        for column_index, column in enumerate(columns):
            if any(
                set(left).intersection(right)
                for left, right in zip(row, column, strict=True)
            ):
                continue
            supports = tuple(
                tuple(sorted(left + right))
                for left, right in zip(row, column, strict=True)
            )
            product_numerator = 1
            for left, right in zip(supports[:-1], supports[1:], strict=True):
                product_numerator *= numerator(plant, left, right)
            matrix[row_index, column_index] = (
                product_numerator / plant.group_size**3
            )
    return matrix


def arbitrary_law_checks(plant: DirectQ4Plant) -> float:
    """Stress every complement/reversal entry under correlated diagonal laws."""

    rng = np.random.default_rng(SEED + 1)
    bound = float((cubic_fixed_pair_energy(Q4) * quintic_fixed_triple_energy(Q4))) ** 0.5
    worst_ratio = 0.0
    for entry in dual_endpoint_entries():
        matrix = occurrence_submatrix(plant, entry, rng)
        laws = [
            (
                np.full(matrix.shape[0], 1 / matrix.shape[0]),
                np.full(matrix.shape[1], 1 / matrix.shape[1]),
            )
        ]
        laws.extend(
            (
                rng.dirichlet(np.full(matrix.shape[0], 0.12)),
                rng.dirichlet(np.full(matrix.shape[1], 0.12)),
            )
            for _ in range(5)
        )
        for row_law, column_law in laws:
            weighted = (
                np.sqrt(row_law)[:, None]
                * matrix
                * np.sqrt(column_law)[None, :]
            )
            nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
            if nuclear > bound * (1 + 2e-12):
                raise AssertionError(("dual endpoint arbitrary law", entry, nuclear, bound))
            worst_ratio = max(worst_ratio, nuclear / bound)
    return worst_ratio


def reversal_checks(plant: DirectQ4Plant) -> int:
    """Check the physical transpose/reversal symmetry, including coordinates."""

    rng = np.random.default_rng(SEED + 2)
    checks = 0
    for left_degree, right_degree in ((3, 1), (1, 5), (5, 3), (3, 3), (1, 1)):
        for _ in range(24):
            left = plant.supports[left_degree][
                int(rng.integers(len(plant.supports[left_degree])))
            ]
            right = plant.supports[right_degree][
                int(rng.integers(len(plant.supports[right_degree])))
            ]
            observed = numerator(plant, left, right)
            reversed_observed = numerator(
                plant,
                transpose_support(right, Q4), transpose_support(left, Q4)
            )
            if observed != reversed_observed:
                raise AssertionError(("physical reversal", left, right))
            checks += 1
    return checks


def selected_q8_link_checks() -> int:
    """Compare direct permutation sums with the exact production evaluator."""

    supports = (
        ((0, 1, 8), (0,)),
        ((0, 1, 2), (9,)),
        ((0,), (0, 1, 2, 3, 4)),
        ((9,), (0, 1, 2, 8, 16)),
        ((6, 9, 19, 41, 50), (7, 12, 46)),
        ((4, 11, 27), (8, 53, 55)),
        ((0,), (9,)),
    )
    nonzero = set()
    for left, right in supports:
        direct = direct_permutation_moment(8, left, right)
        analytic = moment(8, left, right)
        if direct != analytic:
            raise AssertionError(("q8 direct link", left, right, direct, analytic))
        if direct:
            nonzero.add((len(left), len(right)))
    if not {(3, 1), (1, 5), (5, 3), (3, 3), (1, 1)}.issubset(nonzero):
        raise AssertionError(("q8 nonzero link coverage", nonzero))
    return len(supports)


def main() -> None:
    plant = DirectQ4Plant()
    cubic, quintic = exact_q4_endpoint_energies(plant)
    if cubic != Fraction(5, 24):
        raise AssertionError(("q4 cubic endpoint", cubic))
    if quintic != Fraction(7, 8):
        raise AssertionError(("q4 quintic endpoint", quintic))
    if cubic != cubic_fixed_pair_energy(Q4):
        raise AssertionError("cubic formula disagrees with direct enumeration")
    if quintic != quintic_fixed_triple_energy(Q4):
        raise AssertionError("quintic formula disagrees with direct enumeration")

    gram_checks = lifted_gram_checks(plant)
    reversal_count = reversal_checks(plant)
    q8_checks = selected_q8_link_checks()
    worst_ratio = arbitrary_law_checks(plant)
    if len(dual_endpoint_entries()) != 12:
        raise AssertionError("dual endpoint inventory")
    print(
        "q64 dual-endpoint independent audit passed: "
        f"q4_cubic={cubic},q4_quintic={quintic},"
        f"gram_checks={gram_checks},reversal_checks={reversal_count},"
        f"q8_links={q8_checks},arbitrary_law_worst_ratio={worst_ratio:.12g},"
        "verdict=certified"
    )


if __name__ == "__main__":
    main()
