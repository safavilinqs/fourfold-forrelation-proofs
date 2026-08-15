#!/usr/bin/env python3
"""Exact and randomized checks of the global component-assignment dichotomy."""

from __future__ import annotations

from itertools import product

import numpy as np


SEED = 2026071407


def placement_exponent(
    placement: tuple[int, ...],
    components: tuple[tuple[int, ...], ...],
    edge_counts: tuple[int, ...],
    spanning: tuple[bool, ...],
    entries: int,
) -> tuple[str, int]:
    maxima = []
    for component in components:
        counts = [sum(placement[v] == entry for v in component) for entry in range(entries)]
        maxima.append(max(counts))

    if all(maximum <= 1 for maximum in maxima):
        projective_square_exponent = -sum(spanning)
        if projective_square_exponent > -1:
            raise AssertionError(("singleton case lacks spanning suppression", placement, maxima))
        return "projective", projective_square_exponent

    assigned_square_exponent = 0
    for component, edges, maximum in zip(components, edge_counts, maxima):
        exponent = len(component) - maximum - edges
        if exponent > 1 - maximum:
            raise AssertionError(("forest exponent", component, edges, maximum, exponent))
        if exponent > 0:
            raise AssertionError(("weak assignment amplified", component, exponent))
        assigned_square_exponent += exponent
    if assigned_square_exponent > -1:
        raise AssertionError(("strong case lacks suppression", placement, maxima, assigned_square_exponent))
    return "assigned", assigned_square_exponent


def exhaustive_named_graphs() -> tuple[int, int, int]:
    # Two spanning chains, then a spanning chain plus a nonspanning edge.
    examples = (
        (
            (tuple(range(4)), tuple(range(4, 8))),
            (3, 3),
            (True, True),
        ),
        (
            (tuple(range(4)), (4, 5)),
            (3, 1),
            (True, False),
        ),
    )
    total = 0
    projective = 0
    assigned = 0
    for components, edges, spanning in examples:
        vertices = sum(len(component) for component in components)
        for placement in product(range(4), repeat=vertices):
            mode, _ = placement_exponent(placement, components, edges, spanning, 4)
            total += 1
            projective += mode == "projective"
            assigned += mode == "assigned"
    return total, projective, assigned


def randomized_components(trials: int = 100_000) -> None:
    rng = np.random.default_rng(SEED)
    for _ in range(trials):
        component_count = int(rng.integers(1, 6))
        sizes = [int(rng.integers(2, 9)) for _ in range(component_count)]
        starts = np.cumsum([0] + sizes)
        components = tuple(tuple(range(starts[i], starts[i + 1])) for i in range(component_count))
        edges = tuple(size - 1 + int(rng.integers(0, 5)) for size in sizes)
        spanning_flags = [False] * component_count
        spanning_flags[int(rng.integers(0, component_count))] = True
        spanning = tuple(spanning_flags)
        entry_count = int(rng.integers(2, 9))
        placement = tuple(int(value) for value in rng.integers(0, entry_count, size=sum(sizes)))
        placement_exponent(placement, components, edges, spanning, entry_count)


def main() -> None:
    total, projective, assigned = exhaustive_named_graphs()
    randomized_components()
    print(
        "global assignment dichotomy passed: "
        f"exhaustive_placements={total}, projective={projective}, assigned={assigned}, "
        "random_component_instances=100000"
    )


if __name__ == "__main__":
    main()
