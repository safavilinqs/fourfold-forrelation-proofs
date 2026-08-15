#!/usr/bin/env python3
"""Exact dynamic check of the adaptive marked-time dose recurrence."""

from __future__ import annotations

from dataclasses import dataclass
import math
import random


SEED = 2026071414
MAX_MARKS = 12
TOL = 2e-11


@dataclass(frozen=True)
class Node:
    dose: int
    children: tuple["Node", ...] = ()


def branch_dose(node: Node) -> int:
    return node.dose + max((branch_dose(child) for child in node.children), default=0)


def square_root_potential(node: Node) -> float:
    return math.sqrt(node.dose) + max(
        (square_root_potential(child) for child in node.children), default=0.0
    )


def placement_values(node: Node) -> tuple[float, ...]:
    child_values = [placement_values(child) for child in node.children]
    result = [1.0]
    local = 2 * math.sqrt(node.dose)
    for marks in range(1, MAX_MARKS + 1):
        value = 0.0
        for current in range(marks + 1):
            remaining = marks - current
            child = max(
                (values[remaining] for values in child_values),
                default=(1.0 if remaining == 0 else 0.0),
            )
            value += math.comb(marks, current) * local**current * child
        result.append(value)
    return tuple(result)


def random_tree(rng: random.Random, depth: int) -> Node:
    dose = rng.randrange(0, 6)
    if depth == 0:
        return Node(dose)
    width = rng.randrange(0, 7)
    return Node(dose, tuple(random_tree(rng, depth - 1) for _ in range(width)))


def check(node: Node) -> float:
    dose = branch_dose(node)
    potential = square_root_potential(node)
    if potential > dose + TOL:
        raise AssertionError(("sqrt potential exceeds hard dose", potential, dose))
    values = placement_values(node)
    worst = 0.0
    for marks, value in enumerate(values):
        potential_bound = (2 * potential) ** marks if marks else 1.0
        dose_bound = (2 * dose) ** marks if marks else 1.0
        if value > potential_bound * (1 + TOL):
            raise AssertionError(("potential recurrence", marks, value, potential_bound))
        if value > dose_bound * (1 + TOL):
            raise AssertionError(("hard-dose recurrence", marks, value, dose_bound))
        if dose_bound:
            worst = max(worst, value / dose_bound)
    return worst


def main() -> None:
    rng = random.Random(SEED)
    worst = 0.0
    checked = 0
    for _ in range(300):
        node = random_tree(rng, rng.randrange(0, 5))
        worst = max(worst, check(node))
        checked += 1

    # Width duplication must not change a max/complete-frame recurrence.
    child = Node(2, (Node(1), Node(3)))
    narrow = Node(1, (child,))
    wide = Node(1, (child,) * 100)
    if placement_values(narrow) != placement_values(wide):
        raise AssertionError("outcome-width duplication changed the recurrence")

    print(
        "adaptive marked-assignment ledger passed: "
        f"trees={checked}, marks_through={MAX_MARKS}, "
        f"worst_hard_dose_ratio={worst:.12g}, width_duplicate=100"
    )


if __name__ == "__main__":
    main()
