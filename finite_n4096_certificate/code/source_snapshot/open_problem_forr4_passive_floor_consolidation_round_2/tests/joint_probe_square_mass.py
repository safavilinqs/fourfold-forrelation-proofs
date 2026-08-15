#!/usr/bin/env python3
"""Exact checks of the minimal-chain joint-probe square-mass identity."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
import random


SEED = 2026071510
BLOCKS = 4
FULL = (1 << BLOCKS) - 1


def integer_partitions(total: int, maximum: int | None = None):
    if total == 0:
        yield ()
        return
    if maximum is None or maximum > total:
        maximum = total
    for first in range(maximum, 0, -1):
        for rest in integer_partitions(total - first, first):
            yield (first,) + rest


def occupations(dose: int) -> list[tuple[int, ...]]:
    return [
        state
        for state in product(range(dose + 1), repeat=BLOCKS)
        if sum(state) <= dose
    ]


def moment(
    weights: dict[tuple[int, ...], Fraction], mask: int
) -> Fraction:
    result = Fraction(0)
    for state, weight in weights.items():
        value = 1
        for block in range(BLOCKS):
            if (mask >> block) & 1:
                value *= state[block]
        result += weight * value
    return result


def local_polynomial(
    weights: dict[tuple[int, ...], Fraction]
) -> list[Fraction]:
    moments = [moment(weights, mask) for mask in range(1 << BLOCKS)]
    result = [Fraction(0) for _ in range(1 << BLOCKS)]
    for mask in range(1 << BLOCKS):
        subset = mask
        while True:
            result[mask] += moments[subset] * moments[mask ^ subset]
            if subset == 0:
                break
            subset = (subset - 1) & mask
    return result


def local_pair_polynomial(
    weights: dict[tuple[int, ...], Fraction]
) -> list[Fraction]:
    result = [Fraction(0) for _ in range(1 << BLOCKS)]
    for first, first_weight in weights.items():
        for second, second_weight in weights.items():
            weight = first_weight * second_weight
            for mask in range(1 << BLOCKS):
                value = 1
                for block in range(BLOCKS):
                    if (mask >> block) & 1:
                        value *= first[block] + second[block]
                result[mask] += weight * value
    return result


def multiply_square_free(
    left: list[Fraction], right: list[Fraction]
) -> list[Fraction]:
    result = [Fraction(0) for _ in range(1 << BLOCKS)]
    for first, first_value in enumerate(left):
        for second, second_value in enumerate(right):
            if not (first & second):
                result[first | second] += first_value * second_value
    return result


def random_law(
    rng: random.Random, dose: int
) -> dict[tuple[int, ...], Fraction]:
    states = occupations(dose)
    raw = [rng.randrange(1, 1000) for _ in states]
    total = sum(raw)
    return {
        state: Fraction(value, total)
        for state, value in zip(states, raw, strict=True)
    }


def deterministic_identity(
    states: tuple[tuple[int, ...], ...]
) -> Fraction:
    # For deterministic laws S_h=T_h=state_h, (2.1) is pointwise.
    combined = [
        sum(2 * state[block] for state in states)
        for block in range(BLOCKS)
    ]
    total = Fraction(1)
    for value in combined:
        total *= value
    return total


def main() -> None:
    rng = random.Random(SEED)
    checked = 0
    worst_ratio = Fraction(0)
    bound = Fraction(81)

    for dose in range(1, 7):
        law = random_law(rng, dose)
        direct = local_polynomial(law)
        paired = local_pair_polynomial(law)
        if direct != paired:
            raise AssertionError(("local two-copy identity", dose))

    for doses in integer_partitions(6):
        for _ in range(12):
            laws = [random_law(rng, dose) for dose in doses]
            polynomial = [Fraction(0) for _ in range(1 << BLOCKS)]
            polynomial[0] = Fraction(1)
            for law in laws:
                polynomial = multiply_square_free(polynomial, local_polynomial(law))
            value = polynomial[FULL]
            if value > bound:
                raise AssertionError(("joint square-mass bound", doses, value, bound))
            worst_ratio = max(worst_ratio, value / bound)
            checked += 1

    # Pointwise deterministic instances independently check coefficient
    # extraction and the AM-GM input.
    for doses in integer_partitions(6):
        states = tuple(
            occupations(dose)[rng.randrange(len(occupations(dose)))]
            for dose in doses
        )
        value = deterministic_identity(states)
        if value > bound:
            raise AssertionError(("deterministic AM-GM bound", doses, states, value))

    print(
        "joint-probe square mass passed: "
        f"rational_laws={checked}, dose_partitions={sum(1 for _ in integer_partitions(6))}, "
        f"universal_bound={bound}, sqrt_bound=9, "
        f"largest_sample_ratio={float(worst_ratio):.12g}"
    )


if __name__ == "__main__":
    main()
