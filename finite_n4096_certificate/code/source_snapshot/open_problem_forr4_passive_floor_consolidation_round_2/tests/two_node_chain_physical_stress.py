#!/usr/bin/env python3
"""Direct two-node adaptive contractions for every minimal-chain placement.

Each chain vertex is placed in one of root-ket, root-bra, child-ket, or
child-bra.  The root uses an arbitrary complete rank-one frame and each root
outcome selects a new arbitrary complete child frame.  Terminal leaf phases
are chosen adversarially, so the tested value is the sum of absolute leaf
contractions.

This covers a concrete physical marked-sector slice with injective support
maps.  It does not cover unmarked base occupations or larger diagrams.
"""

from __future__ import annotations

from itertools import product

import numpy as np


SEED = 2026071405
TOL = 2e-8


def sylvester(n: int) -> np.ndarray:
    h = np.array([[1.0]])
    while h.shape[0] < n:
        h = np.block([[h, h], [h, -h]])
    return h / np.sqrt(n)


def chain_tensor(n: int) -> np.ndarray:
    h = sylvester(n)
    return np.einsum("ij,jk,kl->ijkl", h, h, h)


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
    if not np.allclose(frame @ frame.conj().T, np.diag(qdiag), atol=2e-12):
        raise AssertionError("frame completeness")
    return frame


def entry_support(vertices: tuple[int, ...], coordinate: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    # The vertex number is also its four-forrelation block, so supports from
    # different vertices cannot collide merely because their numeric labels do.
    return tuple((vertex, coordinate[vertex]) for vertex in vertices)


def node_maps(
    n: int,
    ket_vertices: tuple[int, ...],
    bra_vertices: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray, int]:
    coordinates = tuple(product(range(n), repeat=4))
    supports = set()
    ket_supports = []
    bra_supports = []
    for coordinate in coordinates:
        ket = entry_support(ket_vertices, coordinate)
        bra = entry_support(bra_vertices, coordinate)
        ket_supports.append(ket)
        bra_supports.append(bra)
        supports.add(ket)
        supports.add(bra)
    ordered = sorted(supports, key=lambda value: (len(value), value))
    index = {support: number for number, support in enumerate(ordered)}
    return (
        np.array([index[support] for support in ket_supports], dtype=int),
        np.array([index[support] for support in bra_supports], dtype=int),
        len(ordered),
    )


def one_placement_value(
    rng: np.random.Generator,
    n: int,
    graph: np.ndarray,
    placement: tuple[int, ...],
) -> float:
    entries = [tuple(vertex for vertex, entry in enumerate(placement) if entry == number) for number in range(4)]
    root_ket, root_bra, root_omega = node_maps(n, entries[0], entries[1])
    child_ket, child_bra, child_omega = node_maps(n, entries[2], entries[3])

    root = complete_frame(rng, root_omega)
    graph_flat = graph.reshape(-1)
    lhs = 0.0

    for root_outcome in range(root_omega):
        root_atom = root[root_ket, root_outcome] * np.conj(root[root_bra, root_outcome])
        # This dependence on the root outcome is the adaptive part.
        child = complete_frame(rng, child_omega)
        for child_outcome in range(child_omega):
            child_atom = child[child_ket, child_outcome] * np.conj(
                child[child_bra, child_outcome]
            )
            leaf = np.dot(graph_flat, root_atom * child_atom)
            lhs += abs(leaf)
    return float(lhs)


def stress() -> tuple[float, tuple[int, ...]]:
    rng = np.random.default_rng(SEED)
    n = 2
    graph = chain_tensor(n)
    target = 1 / np.sqrt(n)
    worst = 0.0
    worst_placement = ()

    # Exhaust every placement once, then revisit random placements with new frames.
    placements = list(product(range(4), repeat=4))
    placements.extend(placements[int(rng.integers(0, len(placements)))] for _ in range(512))
    for placement in placements:
        value = one_placement_value(rng, n, graph, placement)
        ratio = value / target
        if ratio > worst:
            worst = ratio
            worst_placement = placement
        if ratio > 1 + TOL:
            raise AssertionError(("two-node chain contraction", ratio, placement))
    return worst, worst_placement


def main() -> None:
    worst, placement = stress()
    print(
        "two-node physical chain stress passed: "
        f"worst_ratio={worst:.12g}, placement={placement}, "
        "exhaustive_placements=256, random_revisits=512"
    )


if __name__ == "__main__":
    main()
