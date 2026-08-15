#!/usr/bin/env python3
"""Direct stress test for assigned/singleton component mixing.

The graph is a spanning four-layer chain plus a disconnected nonspanning
edge.  The chain is forced into the assigned case, while the edge has at
most one vertex in every physical entry.  Arbitrary complete root frames,
root-outcome-dependent child frames, cross-component support entanglement,
and adversarial terminal phases are retained exactly.

Coordinate tuples are restricted to be distinct within each block, as they
are in a Fourier support set.  The test is a finite N=2 falsification slice,
not a proof of the mixed-component induction.
"""

from __future__ import annotations

from itertools import product

import numpy as np


SEED = 2026071406
TOL = 3e-8


LAYERS = (0, 1, 2, 3, 0, 1)
EDGES = ((0, 1), (1, 2), (2, 3), (4, 5))


def sylvester(n: int) -> np.ndarray:
    h = np.array([[1.0]])
    while h.shape[0] < n:
        h = np.block([[h, h], [h, -h]])
    return h / np.sqrt(n)


def valid_coordinates(n: int) -> list[tuple[int, ...]]:
    result = []
    for coordinate in product(range(n), repeat=len(LAYERS)):
        valid = True
        for layer in range(4):
            values = [coordinate[v] for v, block in enumerate(LAYERS) if block == layer]
            if len(values) != len(set(values)):
                valid = False
                break
        if valid:
            result.append(coordinate)
    return result


def graph_values(n: int, coordinates: list[tuple[int, ...]]) -> np.ndarray:
    h = sylvester(n)
    return np.array(
        [np.prod([h[coordinate[left], coordinate[right]] for left, right in EDGES]) for coordinate in coordinates],
        dtype=complex,
    )


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


def support(vertices: tuple[int, ...], coordinate: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    # A physical parity support only sees block and coordinate, not component.
    return tuple(sorted((LAYERS[vertex], coordinate[vertex]) for vertex in vertices))


def node_maps(
    coordinates: list[tuple[int, ...]],
    ket_vertices: tuple[int, ...],
    bra_vertices: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray, int]:
    ket_supports = [support(ket_vertices, coordinate) for coordinate in coordinates]
    bra_supports = [support(bra_vertices, coordinate) for coordinate in coordinates]
    ordered = sorted(set(ket_supports + bra_supports), key=lambda value: (len(value), value))
    index = {value: number for number, value in enumerate(ordered)}
    return (
        np.array([index[value] for value in ket_supports], dtype=int),
        np.array([index[value] for value in bra_supports], dtype=int),
        len(ordered),
    )


def contraction(
    rng: np.random.Generator,
    coordinates: list[tuple[int, ...]],
    graph: np.ndarray,
    placement: tuple[int, ...],
) -> float:
    entries = [tuple(v for v, entry in enumerate(placement) if entry == number) for number in range(4)]
    root_ket, root_bra, root_omega = node_maps(coordinates, entries[0], entries[1])
    child_ket, child_bra, child_omega = node_maps(coordinates, entries[2], entries[3])
    root = complete_frame(rng, root_omega)
    lhs = 0.0

    for root_outcome in range(root_omega):
        root_atom = root[root_ket, root_outcome] * np.conj(root[root_bra, root_outcome])
        child = complete_frame(rng, child_omega)
        for child_outcome in range(child_omega):
            child_atom = child[child_ket, child_outcome] * np.conj(
                child[child_bra, child_outcome]
            )
            lhs += abs(np.dot(graph, root_atom * child_atom))
    return float(lhs)


def admissible(placement: tuple[int, ...]) -> bool:
    chain_counts = [sum(placement[v] == entry for v in range(4)) for entry in range(4)]
    edge_is_singleton = placement[4] != placement[5]
    both_nodes_used = any(entry < 2 for entry in placement) and any(entry >= 2 for entry in placement)
    return max(chain_counts) >= 2 and edge_is_singleton and both_nodes_used


def stress() -> tuple[float, tuple[int, ...], int]:
    rng = np.random.default_rng(SEED)
    n = 2
    coordinates = valid_coordinates(n)
    graph = graph_values(n, coordinates)
    target = 1 / np.sqrt(n)
    placements = [placement for placement in product(range(4), repeat=6) if admissible(placement)]

    # Cover every admissible placement once.  New frames are drawn per placement.
    worst = 0.0
    worst_placement = ()
    for placement in placements:
        value = contraction(rng, coordinates, graph, placement)
        ratio = value / target
        if ratio > worst:
            worst = ratio
            worst_placement = placement
        if ratio > 1 + TOL:
            raise AssertionError(("mixed assigned/singleton contraction", ratio, placement))
    return worst, worst_placement, len(placements)


def main() -> None:
    worst, placement, checked = stress()
    print(
        "mixed-component physical stress passed: "
        f"worst_ratio={worst:.12g}, placement={placement}, "
        f"placements={checked}, distinct_coordinate_tuples=4"
    )


if __name__ == "__main__":
    main()
