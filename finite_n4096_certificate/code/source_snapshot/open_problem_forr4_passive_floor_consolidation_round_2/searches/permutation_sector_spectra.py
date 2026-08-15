#!/usr/bin/env python3
"""Spectra of pure signed-permutation match sectors.

For r distinct odd column labels on the X side and r distinct odd row labels
on the Y side, the exact moment matrix is q^r/(q)_r times the permanent
compound of the single-particle H_{q^2}.  This script verifies the identity
and reports small-q singular spectra.
"""

from __future__ import annotations

from itertools import combinations, permutations, product
from math import factorial

import numpy as np


def sylvester(q: int) -> np.ndarray:
    h = np.array([[1.0]])
    while h.shape[0] < q:
        h = np.block([[h, h], [h, -h]])
    return h / np.sqrt(q)


def falling(q: int, r: int) -> int:
    value = 1
    for offset in range(r):
        value *= q - offset
    return value


def sector_states(q: int, side: str, r: int):
    states = []
    for labels in combinations(range(q), r):
        for characters in product(range(q), repeat=r):
            if side == "left":
                states.append(tuple((characters[i], labels[i]) for i in range(r)))
            else:
                states.append(tuple((labels[i], characters[i]) for i in range(r)))
    return states


def permanent(matrix: np.ndarray) -> float:
    r = matrix.shape[0]
    return float(sum(np.prod([matrix[i, perm[i]] for i in range(r)]) for perm in permutations(range(r))))


def exact_moment(left, right, k):
    r = len(left)
    total = 0.0
    for perm in permutations(range(r)):
        phase = 1.0
        for i in range(r):
            x, y = left[i]
            u, v = right[perm[i]]
            phase *= (np.sqrt(k.shape[0]) * k[x, u]) * (np.sqrt(k.shape[0]) * k[v, y])
        total += phase
    return total / falling(k.shape[0], r)


def sector_matrix(q: int, r: int):
    k = sylvester(q)
    h = np.kron(k, k)
    left_states = sector_states(q, "left", r)
    right_states = sector_states(q, "right", r)
    exact = np.empty((len(left_states), len(right_states)))
    compound = np.empty_like(exact)
    scale = q**r / falling(q, r)

    for i, left in enumerate(left_states):
        left_indices = [x * q + y for x, y in left]
        for j, right in enumerate(right_states):
            right_indices = [u * q + v for u, v in right]
            exact[i, j] = exact_moment(left, right, k)
            compound[i, j] = scale * permanent(h[np.ix_(left_indices, right_indices)])
    residual = float(np.max(np.abs(exact - compound)))
    return exact, residual


def main() -> None:
    rows = []
    for q in (2, 4):
        for r in range(1, q + 1):
            matrix, residual = sector_matrix(q, r)
            singular = np.linalg.svd(matrix, compute_uv=False)
            rank = int(np.sum(singular > 2e-10))
            op = float(singular[0])
            fro = float(np.linalg.norm(singular))
            upper = q**r / falling(q, r)
            if residual > 2e-12:
                raise AssertionError(("bosonic compound identity", q, r, residual))
            if op > upper * (1 + 2e-10):
                raise AssertionError(("compression norm", q, r, op, upper))
            rows.append(
                f"q={q},r={r},dim={matrix.shape[0]},rank={rank},"
                f"op={op:.12g},fro={fro:.12g},compound_cap={upper:.12g}"
            )
    print("signed-permutation sector spectra:\n" + "\n".join(rows))


if __name__ == "__main__":
    main()
