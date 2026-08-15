#!/usr/bin/env python3
"""N=4 spot checks for the strong-plus-weak assignment side of the dichotomy."""

from __future__ import annotations

from itertools import product

import numpy as np


SEED = 2026071409
TOL = 4e-8
LAYERS = (0, 1, 2, 3, 0, 1)
EDGES = ((0, 1), (1, 2), (2, 3), (4, 5))


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


def complete_frame(rng, omega):
    qdiag = rng.random(omega)
    qdiag /= qdiag.sum()
    return np.sqrt(qdiag)[:, None] * random_unitary(rng, omega)


def valid_coordinates(n):
    return [
        coordinate
        for coordinate in product(range(n), repeat=6)
        if coordinate[0] != coordinate[4] and coordinate[1] != coordinate[5]
    ]


def graph_values(n, coordinates):
    h = sylvester(n)
    return np.array(
        [np.prod([h[c[left], c[right]] for left, right in EDGES]) for c in coordinates],
        dtype=complex,
    )


def support(vertices, coordinate):
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
    n = 4
    coordinates = valid_coordinates(n)
    graph = graph_values(n, coordinates)
    target = 1 / np.sqrt(n)
    placements = (
        (0, 0, 2, 3, 1, 2),
        (1, 1, 2, 3, 0, 3),
        (0, 1, 2, 2, 3, 0),
        (0, 1, 3, 3, 2, 1),
        (2, 2, 0, 1, 3, 0),
        (3, 3, 0, 1, 2, 1),
    )
    worst = 0.0
    worst_placement = ()
    for placement in placements:
        value = contraction(rng, coordinates, graph, placement)
        ratio = value / target
        if ratio > worst:
            worst = ratio
            worst_placement = placement
        if ratio > 1 + TOL:
            raise AssertionError(("N=4 mixed contraction", ratio, placement))
    print(
        "N=4 mixed-component spot checks passed: "
        f"worst_ratio={worst:.12g}, placement={worst_placement}, "
        f"placements={len(placements)}, distinct_coordinate_tuples={len(coordinates)}"
    )


if __name__ == "__main__":
    main()
