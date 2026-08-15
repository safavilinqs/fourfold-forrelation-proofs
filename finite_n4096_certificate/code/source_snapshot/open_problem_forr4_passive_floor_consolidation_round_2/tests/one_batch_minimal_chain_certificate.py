#!/usr/bin/env python3
"""Exact KKT certificate for the dose-six one-batch minimal-chain bound."""

from __future__ import annotations

from itertools import product

import sympy as sp


N = sp.Integer(1024)
DOSE = 6
EDGES = {(0, 1), (1, 2), (2, 3)}
OPTIMAL_SUPPORT = (
    (1, 2, 2, 1),
    (2, 1, 1, 2),
    (1, 2, 1, 2),
    (2, 1, 2, 1),
)


def feature(state: tuple[int, ...], mask: int) -> sp.Integer:
    value = sp.Integer(1)
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
    rows = [
        [int((min(u, v), max(u, v)) in EDGES) for v in right]
        for u in left
    ]
    return gf2_rank(rows)


def main() -> None:
    states = [state for state in product(range(DOSE + 1), repeat=4) if sum(state) <= DOSE]
    moments = [
        sum(feature(state, mask) for state in OPTIMAL_SUPPORT) / sp.Integer(4)
        for mask in range(16)
    ]
    weights = [N ** sp.Rational(cut_rank(mask) - 2, 2) for mask in range(16)]
    objective = sp.simplify(
        sum(
            weights[mask] * sp.sqrt(moments[mask] * moments[15 ^ mask])
            for mask in range(16)
        )
    )
    expected = sp.Rational(2337, 256) + 3 * sp.sqrt(2) / 8
    if sp.simplify(objective - expected) != 0:
        raise AssertionError(("objective", objective, expected))

    equality_states = []
    minimum_gap = None
    next_state = None
    for state in states:
        gradient = sp.Integer(0)
        for mask in range(16):
            complement = 15 ^ mask
            gradient += weights[mask] * sp.Rational(1, 2) * (
                sp.sqrt(moments[complement] / moments[mask]) * feature(state, mask)
                + sp.sqrt(moments[mask] / moments[complement]) * feature(state, complement)
            )
        gap = sp.simplify(objective - gradient)
        if gap == 0:
            equality_states.append(state)
        elif gap.is_positive is not True:
            raise AssertionError(("KKT gradient violation or undecided sign", state, gradient, gap))
        if gap != 0 and (minimum_gap is None or float(gap) < float(minimum_gap)):
            minimum_gap = gap
            next_state = state

    if tuple(equality_states) != tuple(sorted(OPTIMAL_SUPPORT)):
        raise AssertionError(("KKT equality support", equality_states, OPTIMAL_SUPPORT))
    expected_gap = sp.Rational(33, 512) + 3 * sp.sqrt(2) / 32
    if sp.simplify(minimum_gap - expected_gap) != 0:
        raise AssertionError(("next KKT gap", minimum_gap, expected_gap, next_state))

    threshold_gap = sp.simplify(sp.Rational(32, 3) - objective)
    if threshold_gap.is_positive is not True:
        raise AssertionError(("N=1024 threshold", objective, threshold_gap))

    print(
        "one-batch minimal-chain certificate passed: "
        f"states={len(states)}, objective={objective}, value={float(objective):.12g}, "
        f"TV_bound={float(objective / 32):.12g}, "
        f"threshold_margin={float(sp.Rational(1, 3) - objective / 32):.12g}, "
        f"next_KKT_gap={minimum_gap}, next_state={next_state}"
    )


if __name__ == "__main__":
    main()
