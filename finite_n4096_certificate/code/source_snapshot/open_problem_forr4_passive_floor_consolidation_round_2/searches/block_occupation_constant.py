#!/usr/bin/env python3
"""Optimize the exact block-occupation frame constant for the minimal chain.

For a probe distribution over four block occupations n with sum at most t,
the objective sums the ket/bra Cauchy masses over all sixteen assignments of
the four distinct block marks.  The raw objective charges every cut equally.
The spectral objective weights each cut by its exact Hadamard-chain rank
factor relative to N^{-1/2}.  Both are concave and solved with
CVXPY/CLARABEL.  Values are diagnostics until paired with rigorous dual
certificates.
"""

from __future__ import annotations

from itertools import product

import cvxpy as cp
import numpy as np


def occupations(dose: int) -> list[tuple[int, int, int, int]]:
    return [state for state in product(range(dose + 1), repeat=4) if sum(state) <= dose]


def feature(state: tuple[int, ...], mask: int) -> int:
    value = 1
    for block in range(4):
        if (mask >> block) & 1:
            value *= state[block]
    return value


def gf2_rank(rows: list[list[int]]) -> int:
    if not rows:
        return 0
    width = len(rows[0])
    values = [sum((bit & 1) << column for column, bit in enumerate(row)) for row in rows]
    rank = 0
    for column in range(width):
        pivot = next((i for i in range(rank, len(values)) if (values[i] >> column) & 1), None)
        if pivot is None:
            continue
        values[rank], values[pivot] = values[pivot], values[rank]
        for i in range(len(values)):
            if i != rank and ((values[i] >> column) & 1):
                values[i] ^= values[rank]
        rank += 1
    return rank


def cut_rank(mask: int) -> int:
    left = [v for v in range(4) if (mask >> v) & 1]
    right = [v for v in range(4) if not ((mask >> v) & 1)]
    edges = {(0, 1), (1, 2), (2, 3)}
    rows = []
    for u in left:
        row = []
        for v in right:
            row.append(int((min(u, v), max(u, v)) in edges))
        rows.append(row)
    return gf2_rank(rows)


def solve(
    dose: int, n_dimension: int | None = None
) -> tuple[float, np.ndarray, list[tuple[int, int, int, int]]]:
    states = occupations(dose)
    features = np.array([[feature(state, mask) for state in states] for mask in range(16)], dtype=float)
    weights = cp.Variable(len(states), nonneg=True)
    moments = features @ weights
    terms = []
    # Each mask and its complement are both physical ket/bra orientations.
    for mask in range(16):
        complement = 15 ^ mask
        weight = 1.0
        if n_dimension is not None:
            # Weighted nuclear scale N^{(r-e)/2}, e=3.  Divide by N^{-1/2}
            # so the objective is the explicit coefficient of 1/sqrt(N).
            rank = cut_rank(mask)
            weight = n_dimension ** ((rank - 2) / 2)
        terms.append(weight * cp.geo_mean(cp.hstack([moments[mask], moments[complement]])))
    problem = cp.Problem(cp.Maximize(cp.sum(terms)), [cp.sum(weights) == 1])
    value = problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-10,
        tol_gap_rel=1e-10,
        tol_feas=1e-10,
        max_iter=1000,
    )
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise AssertionError(("occupation optimization", dose, problem.status))
    return float(value), np.asarray(weights.value), states


def main() -> None:
    rows = []
    weighted_dose_six = None
    for dose in range(2, 7):
        raw_value, _, _ = solve(dose)
        value, weights, states = solve(dose, n_dimension=1024)
        support = sorted(
            ((float(weight), state) for weight, state in zip(weights, states) if weight > 2e-7),
            reverse=True,
        )
        rows.append(
            f"dose={dose},raw={raw_value:.12g},N1024_spectral={value:.12g},support="
            + ";".join(f"{state}:{weight:.6g}" for weight, state in support[:12])
        )
        if dose == 2 and not np.isclose(value, 1.0, atol=3e-6):
            # The optimizer puts its sector mass on a crossing pairing.
            raise AssertionError(("dose-two spectral normalization", value))
        if dose == 6:
            weighted_dose_six = value
    if weighted_dose_six is None or weighted_dose_six >= 32 / 3:
        raise AssertionError(
            ("dose-six one-batch spectral constant misses N=1024 threshold", weighted_dose_six)
        )
    print("minimal-chain block-occupation constants:\n" + "\n".join(rows))
    print(f"dose-six spectral threshold slack={(32/3)/weighted_dose_six:.12g}")


if __name__ == "__main__":
    main()
