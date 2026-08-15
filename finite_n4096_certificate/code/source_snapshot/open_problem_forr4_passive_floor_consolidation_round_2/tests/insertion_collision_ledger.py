#!/usr/bin/env python3
"""Exact enumeration of the combined parity-collision/Bessel ledger."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, permutations
import math
import random


SEED = 2026071413


def falling(value: int, order: int) -> int:
    result = 1
    for offset in range(order):
        result *= value - offset
    return result


def supports(universe: tuple[int, ...], dose: int) -> list[frozenset[int]]:
    return [
        frozenset(chosen)
        for size in range(dose + 1)
        for chosen in combinations(universe, size)
    ]


def one_entry_check(
    universe: tuple[int, ...], dose: int, marks: int, weights: dict[frozenset[int], Fraction]
) -> tuple[Fraction, int]:
    direct = Fraction(0)
    fibers: dict[tuple[frozenset[int], frozenset[int]], int] = {}

    for support, weight in weights.items():
        for selected_set in combinations(sorted(support), marks):
            base = support.difference(selected_set)
            for ordered in permutations(selected_set):
                recovered = base.union(ordered)
                if recovered != support:
                    raise AssertionError(("support reconstruction", support, base, ordered))
                direct += weight
                key = (frozenset(base), frozenset(recovered))
                fibers[key] = fibers.get(key, 0) + 1

    formula = sum(
        (weight * falling(len(support), marks) for support, weight in weights.items()),
        Fraction(0),
    )
    if direct != formula:
        raise AssertionError(("falling-factorial identity", dose, marks, direct, formula))
    if formula > dose**marks:
        raise AssertionError(("dose bound", dose, marks, formula, dose**marks))

    maximum_fiber = max(fibers.values(), default=1)
    if maximum_fiber > math.factorial(marks):
        raise AssertionError(("collision fiber", dose, marks, maximum_fiber))
    if marks and any(len(support) >= marks for support in weights):
        if maximum_fiber != math.factorial(marks):
            raise AssertionError(("factorial not attained", dose, marks, maximum_fiber))
    return formula, maximum_fiber


def main() -> None:
    rng = random.Random(SEED)
    universe = tuple(range(7))
    checked = 0
    worst_fraction = Fraction(0)
    largest_fiber = 0

    for dose in range(0, 7):
        family = supports(universe, dose)
        raw = [rng.randrange(1, 1000) for _ in family]
        total = sum(raw)
        weights = {
            support: Fraction(weight, total)
            for support, weight in zip(family, raw, strict=True)
        }
        for marks in range(0, min(4, dose) + 1):
            value, fiber = one_entry_check(universe, dose, marks, weights)
            scale = Fraction(dose**marks) if dose or marks == 0 else Fraction(1)
            fraction = value / scale if scale else Fraction(0)
            worst_fraction = max(worst_fraction, fraction)
            largest_fiber = max(largest_fiber, fiber)
            checked += 1

    print(
        "combined insertion/collision ledger passed: "
        f"exact_cases={checked}, largest_fiber={largest_fiber}, "
        f"worst_dose_fraction={float(worst_fraction):.12g}"
    )


if __name__ == "__main__":
    main()
