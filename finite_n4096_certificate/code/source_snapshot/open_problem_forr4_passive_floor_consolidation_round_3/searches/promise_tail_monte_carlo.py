#!/usr/bin/env python3
"""Seeded promise-tail diagnostic for attenuated exact plants.

This Monte Carlo calculation ranks concentration as a repair route after the
mixed-orbit obstruction.  It is explicitly not a rigorous tail bound.
"""

from __future__ import annotations

import argparse

import numpy as np


SEED = 20260716


def sylvester_sign(order: int) -> np.ndarray:
    result = np.asarray([[1]], dtype=np.int8)
    while len(result) < order:
        result = np.block([[result, result], [result, -result]])
    return result


def signed_permutation_pair(
    hadamard: np.ndarray, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray]:
    order = len(hadamard)
    permutation = rng.permutation(order)
    signs = rng.choice((-1, 1), order)
    left = (hadamard[:, permutation] * signs).reshape(-1)
    right = np.empty((order, order), dtype=np.int8)
    right[permutation, :] = signs[:, None] * hadamard
    return left, right.reshape(-1)


def exact_plant(
    hadamard: np.ndarray, rng: np.random.Generator
) -> tuple[np.ndarray, ...]:
    first = signed_permutation_pair(hadamard, rng)
    second = signed_permutation_pair(hadamard, rng)
    third = signed_permutation_pair(hadamard, rng)
    return (
        first[0],
        first[1] * second[0],
        second[1] * third[0],
        third[1],
    )


def normalized_fwht(values: np.ndarray) -> np.ndarray:
    result = np.asarray(values, dtype=float).copy()
    dimension = result.shape[1]
    width = 1
    while width < dimension:
        blocks = result.reshape(len(result), -1, 2 * width)
        left = blocks[:, :, :width].copy()
        right = blocks[:, :, width:].copy()
        blocks[:, :, :width] = left + right
        blocks[:, :, width:] = left - right
        width *= 2
    return result / np.sqrt(dimension)


def four_chain(blocks: tuple[np.ndarray, ...]) -> np.ndarray:
    dimension = blocks[0].shape[1]
    state = normalized_fwht(blocks[3])
    state *= blocks[2]
    state = normalized_fwht(state)
    state *= blocks[1]
    state = normalized_fwht(state)
    return np.sum(blocks[0] * state, axis=1) / dimension


def sample_batch(
    plant: tuple[np.ndarray, ...],
    beta: float,
    batch: int,
    rng: np.random.Generator,
) -> np.ndarray:
    positive_probability = (1 + beta) / 2
    blocks = tuple(
        values
        * (2 * (rng.random((batch, len(values))) < positive_probability) - 1)
        for values in plant
    )
    return four_chain(blocks)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=32)
    parser.add_argument("--beta", type=float, default=0.780899845855353)
    parser.add_argument("--plants", type=int, default=25)
    parser.add_argument("--samples-per-plant", type=int, default=4000)
    parser.add_argument("--batch", type=int, default=400)
    arguments = parser.parse_args()
    if arguments.samples_per_plant % arguments.batch:
        raise ValueError("samples per plant must be divisible by batch")

    rng = np.random.default_rng(SEED)
    hadamard = sylvester_sign(arguments.order)
    failures = 0
    total = 0
    minimum = 1.0
    for _ in range(arguments.plants):
        plant = exact_plant(hadamard, rng)
        exact = four_chain(
            tuple(values[None, :] for values in plant)
        )[0]
        if not np.isclose(exact, 1, atol=2e-13):
            raise AssertionError(("exact plant", exact))
        for _ in range(arguments.samples_per_plant // arguments.batch):
            values = sample_batch(
                plant, arguments.beta, arguments.batch, rng
            )
            failures += int(np.count_nonzero(values < 1 / 4))
            total += len(values)
            minimum = min(minimum, float(values.min()))
    print(
        "promise-tail Monte Carlo diagnostic: "
        f"seed={SEED},order={arguments.order},N={arguments.order**2},"
        f"beta={arguments.beta:.15g},samples={total},"
        f"failures={failures},rate={failures/total:.12g},"
        f"minimum={minimum:.12g},status=NONRIGOROUS"
    )


if __name__ == "__main__":
    main()
