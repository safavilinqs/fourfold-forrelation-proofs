#!/usr/bin/env python3
"""Physical stress test for the all-singleton/projective side of the dichotomy."""

from __future__ import annotations

from itertools import permutations, product

import numpy as np


SEED = 2026071408
TOL = 3e-8
LAYERS = (0, 1, 2, 3, 0, 1, 2, 3)
EDGES = ((0, 1), (1, 2), (2, 3), (4, 5), (5, 6), (6, 7))


def sylvester(n: int) -> np.ndarray:
    h = np.array([[1.0]])
    while h.shape[0] < n:
        h = np.block([[h, h], [h, -h]])
    return h / np.sqrt(n)


def random_unitary(rng: np.random.Generator, n: int) -> np.ndarray:
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, r = np.linalg.qr(z)
    diagonal = np.diag(r)
    phase = np.where(np.abs(diagonal) > 0, diagonal / np.abs(diagonal), 1.0)
    return q * phase.conj()


def complete_frame(rng: np.random.Generator, omega: int) -> np.ndarray:
    qdiag = rng.random(omega)
    qdiag /= qdiag.sum()
    frame = np.sqrt(qdiag)[:, None] * random_unitary(rng, omega)
    return frame


def valid_coordinates(n: int) -> list[tuple[int, ...]]:
    result = []
    for coordinate in product(range(n), repeat=8):
        if all(
            len({coordinate[v] for v, block in enumerate(LAYERS) if block == layer}) == 2
            for layer in range(4)
        ):
            result.append(coordinate)
    return result


def graph_values(n: int, coordinates: list[tuple[int, ...]]) -> np.ndarray:
    h = sylvester(n)
    return np.array(
        [np.prod([h[c[left], c[right]] for left, right in EDGES]) for c in coordinates],
        dtype=complex,
    )


def support(vertices: tuple[int, ...], coordinate: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted((LAYERS[v], coordinate[v]) for v in vertices))


def node_maps(coordinates, ket_vertices, bra_vertices):
    ket = [support(ket_vertices, c) for c in coordinates]
    bra = [support(bra_vertices, c) for c in coordinates]
    universe = sorted(set(ket + bra), key=lambda value: (len(value), value))
    index = {value: number for number, value in enumerate(universe)}
    return np.array([index[v] for v in ket]), np.array([index[v] for v in bra]), len(universe)


def contraction(rng, coordinates, graph, placement):
    entries = [tuple(v for v, entry in enumerate(placement) if entry == number) for number in range(4)]
    rk, rb, ro = node_maps(coordinates, entries[0], entries[1])
    ck, cb, co = node_maps(coordinates, entries[2], entries[3])
    root = complete_frame(rng, ro)
    lhs = 0.0
    for y in range(ro):
        root_atom = root[rk, y] * np.conj(root[rb, y])
        child = complete_frame(rng, co)
        for z in range(co):
            child_atom = child[ck, z] * np.conj(child[cb, z])
            lhs += abs(np.dot(graph, root_atom * child_atom))
    return float(lhs)


def main() -> None:
    rng = np.random.default_rng(SEED)
    n = 2
    coordinates = valid_coordinates(n)
    graph = graph_values(n, coordinates)
    target = 1 / n  # two spanning singleton-projective components
    worst = 0.0
    worst_placement = ()
    checked = 0
    first_component = (0, 1, 2, 3)
    for second_entries in permutations(range(4)):
        placement = first_component + second_entries
        for _ in range(4):
            value = contraction(rng, coordinates, graph, placement)
            ratio = value / target
            if ratio > worst:
                worst = ratio
                worst_placement = placement
            if ratio > 1 + TOL:
                raise AssertionError(("all-singleton contraction", ratio, placement))
            checked += 1
    print(
        "all-singleton physical stress passed: "
        f"worst_ratio={worst:.12g}, placement={worst_placement}, "
        f"permutation_frame_instances={checked}, distinct_coordinate_tuples={len(coordinates)}"
    )


if __name__ == "__main__":
    main()
