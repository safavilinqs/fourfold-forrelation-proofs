#!/usr/bin/env python3
"""Derivative-energy and finite-size checks for the attenuated exact plant."""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations, product

import numpy as np


SEED = 2026071411


def sylvester_sign(order: int) -> np.ndarray:
    result = np.array([[1]], dtype=int)
    while result.shape[0] < order:
        result = np.block([[result, result], [result, -result]])
    if result.shape != (order, order):
        raise ValueError("order must be a power of two")
    return result


def signed_permutation(
    order: int, permutation: tuple[int, ...], signs: tuple[int, ...]
) -> np.ndarray:
    result = np.zeros((order, order), dtype=int)
    for column, row in enumerate(permutation):
        result[row, column] = signs[column]
    return result


def planted_pair(k: np.ndarray, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return (k @ p).reshape(-1), (p @ k).reshape(-1)


def plus_plant(
    pairs: tuple[tuple[np.ndarray, np.ndarray], ...]
) -> tuple[np.ndarray, ...]:
    (x1, y1), (x2, y2), (x3, y3) = pairs
    return x1, y1 * x2, y2 * x3, y3


def chain_tensor(hadamard: np.ndarray) -> np.ndarray:
    dimension = hadamard.shape[0]
    return np.einsum(
        "ij,jk,kl->ijkl", hadamard, hadamard, hadamard, optimize=True
    ) / dimension


def derivative_energies(
    tensor: np.ndarray, blocks: tuple[np.ndarray, ...]
) -> tuple[float, dict[int, float]]:
    dimension = len(blocks[0])
    dressed = tensor.copy()
    for block, values in enumerate(blocks):
        shape = [1] * 4
        shape[block] = dimension
        dressed *= values.reshape(shape)

    value = float(dressed.sum())
    energies: dict[int, float] = {}
    for mask in range(1, 1 << 4):
        open_blocks = {block for block in range(4) if (mask >> block) & 1}
        derivative = dressed
        for axis in sorted(set(range(4)) - open_blocks, reverse=True):
            derivative = derivative.sum(axis=axis)
        energies[mask] = float(np.sum(derivative * derivative))
    return value, energies


def check_plant(
    tensor: np.ndarray, blocks: tuple[np.ndarray, ...], tolerance: float = 1e-12
) -> float:
    dimension = len(blocks[0])
    value, energies = derivative_energies(tensor, blocks)
    if abs(value - 1) > tolerance:
        raise AssertionError(("planted value", dimension, value))
    target = 1 / dimension
    error = max(abs(energy - target) for energy in energies.values())
    if error > tolerance:
        raise AssertionError(("derivative energy", dimension, error, energies))
    return error


def exhaustive_q2() -> tuple[int, float]:
    order = 2
    k = sylvester_sign(order)
    hadamard = np.kron(k, k) / order
    tensor = chain_tensor(hadamard)
    pairs = [
        planted_pair(k, signed_permutation(order, permutation, signs))
        for permutation in permutations(range(order))
        for signs in product((-1, 1), repeat=order)
    ]
    count = 0
    maximum_error = 0.0
    for triple in product(pairs, repeat=3):
        maximum_error = max(maximum_error, check_plant(tensor, plus_plant(triple)))
        count += 1
    return count, maximum_error


def random_q4(trials: int = 100) -> float:
    rng = np.random.default_rng(SEED)
    order = 4
    k = sylvester_sign(order)
    hadamard = np.kron(k, k) / order
    tensor = chain_tensor(hadamard)
    maximum_error = 0.0
    for _ in range(trials):
        pairs = []
        for _ in range(3):
            permutation = tuple(int(x) for x in rng.permutation(order))
            signs = tuple(int(x) for x in rng.choice((-1, 1), order))
            pairs.append(
                planted_pair(k, signed_permutation(order, permutation, signs))
            )
        maximum_error = max(
            maximum_error, check_plant(tensor, plus_plant(tuple(pairs)))
        )
    return maximum_error


def finite_size_budget() -> dict[str, Fraction]:
    beta = Fraction(5, 6)
    dimension = 1024
    mean = beta**4
    variance = (1 - beta**8) / dimension
    promise_gap = mean - Fraction(1, 4)
    cantelli_failure = variance / (variance + promise_gap**2)
    minimal_sector = Fraction(9, 32) * beta**4
    rejected_hypothetical_budget = minimal_sector + 2 * cantelli_failure
    if rejected_hypothetical_budget >= Fraction(1, 3):
        raise AssertionError(
            ("hypothetical arithmetic changed", rejected_hypothetical_budget)
        )
    return {
        "beta": beta,
        "mean": mean,
        "variance": variance,
        "cantelli_failure": cantelli_failure,
        "minimal_sector": minimal_sector,
        "rejected_hypothetical_budget": rejected_hypothetical_budget,
    }


def main() -> None:
    exhaustive_count, q2_error = exhaustive_q2()
    q4_error = random_q4()
    budget = finite_size_budget()
    print(
        "attenuated exact plant passed: "
        f"q2_plants={exhaustive_count}, q2_max_error={q2_error:.3g}, "
        f"q4_random_plants=100, q4_max_error={q4_error:.3g}, "
        f"beta={budget['beta']}, conditional_variance={float(budget['variance']):.12g}, "
        f"cantelli_failure={float(budget['cantelli_failure']):.12g}, "
        "near_unit_square_function_budget=REJECTED"
    )


if __name__ == "__main__":
    main()
