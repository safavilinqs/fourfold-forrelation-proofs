#!/usr/bin/env python3
"""Dimension-sensitive trace norm for two nonadaptive dose-one probes."""

from __future__ import annotations

from itertools import permutations

import numpy as np


SEED = 2026071417
SQUARE_MASS = 3 / 32


def sylvester(order: int) -> np.ndarray:
    result = np.array([[1.0]])
    while result.shape[0] < order:
        result = np.block([[result, result], [result, -result]])
    return result / np.sqrt(order)


def uniform_minimal_operator(order: int) -> np.ndarray:
    hadamard = sylvester(order)
    chain = np.einsum(
        "ij,jk,kl->ijkl", hadamard, hadamard, hadamard, optimize=True
    )
    modes = 4 * order
    result = np.zeros((modes * modes, modes * modes))
    weight = 1 / modes**2
    for block_order in permutations(range(4)):
        first, second, third, fourth = block_order
        rows = np.array(
            [
                (first * order + i) * modes + (second * order + j)
                for i in range(order)
                for j in range(order)
            ]
        )
        columns = np.array(
            [
                (third * order + k) * modes + (fourth * order + ell)
                for k in range(order)
                for ell in range(order)
            ]
        )
        flattening = np.transpose(chain, axes=block_order).reshape(
            order**2, order**2
        )
        result[np.ix_(rows, columns)] = weight * flattening
    return result


def analytic_trace_norm(order: int) -> float:
    return 1 / (2 * np.sqrt(order)) + 1 / (4 * order)


def block_mass_bound(
    order: int, first: np.ndarray, second: np.ndarray
) -> float:
    crossing = 1 / np.sqrt(order)
    adjacent = 1 / order
    partitions = (
        ((0, 2), (1, 3), crossing),
        ((0, 3), (1, 2), crossing),
        ((0, 1), (2, 3), adjacent),
    )
    result = 0.0
    for left, right, weight in partitions:
        result += 2 * weight * (
            np.sqrt(
                first[left[0]]
                * first[left[1]]
                * second[right[0]]
                * second[right[1]]
            )
            + np.sqrt(
                first[right[0]]
                * first[right[1]]
                * second[left[0]]
                * second[left[1]]
            )
        )
    return float(result)


def main() -> None:
    rows = []
    for order in (1, 2, 4, 8):
        operator = uniform_minimal_operator(order)
        if not np.allclose(operator, operator.T, atol=2e-13):
            raise AssertionError(("Hermiticity", order))
        trace_norm = float(np.abs(np.linalg.eigvalsh(operator)).sum())
        expected = analytic_trace_norm(order)
        if not np.isclose(trace_norm, expected, atol=2e-12):
            raise AssertionError(("trace norm formula", order, trace_norm, expected))
        square_ratio = trace_norm * np.sqrt(order / SQUARE_MASS)
        rows.append(f"N={order}:trace={trace_norm:.12g},ratio={square_ratio:.12g}")

    rng = np.random.default_rng(SEED)
    checked = 0
    for order in (2, 4, 16, 1024):
        upper = analytic_trace_norm(order)
        uniform = np.full(4, 1 / 4)
        if not np.isclose(block_mass_bound(order, uniform, uniform), upper):
            raise AssertionError(("uniform equality", order))
        for _ in range(20_000):
            first = rng.dirichlet(np.ones(4) * 10 ** rng.uniform(-1, 1))
            second = rng.dirichlet(np.ones(4) * 10 ** rng.uniform(-1, 1))
            value = block_mass_bound(order, first, second)
            if value > upper * (1 + 2e-12):
                raise AssertionError(("block-mass optimization", order, value, upper))
            checked += 1

    limit_ratio = 0.5 / np.sqrt(SQUARE_MASS)
    if not np.isclose(limit_ratio, np.sqrt(8 / 3)):
        raise AssertionError(("limiting ratio", limit_ratio))
    print(
        "two-batch trace norm passed: "
        + ", ".join(rows)
        + f", random_mass_pairs={checked}, limiting_square_constant={limit_ratio:.12g}"
    )


if __name__ == "__main__":
    main()
