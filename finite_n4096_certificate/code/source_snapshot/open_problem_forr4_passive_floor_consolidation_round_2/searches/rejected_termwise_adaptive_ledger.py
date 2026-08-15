#!/usr/bin/env python3
"""Reject the termwise-l1 adaptive minimal-chain ledger at hard dose six.

For each node dose t and local ket/bra mark counts (a,b), a concave program
maximizes the exact block-occupation Bessel mass.  Every placement of the
four distinct block marks into the ket/bra entries of every integer dose
partition is then summed.  Each placement receives the repaired graph factor
N^{(1-k_max)/2}, or N^{-1/2} in the all-singleton case.

The local probe is optimized separately for every placement, so this is a
safe but very loose ledger bound.  The script records that it cannot support
the realistic-size theorem and identifies the all-one-photon obstruction.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import product

import cvxpy as cp
import numpy as np


N = 1024
DOSE = 6


def occupations(t: int) -> list[tuple[int, int, int, int]]:
    return [state for state in product(range(t + 1), repeat=4) if sum(state) <= t]


def monomial(state: tuple[int, ...], blocks: range) -> int:
    value = 1
    for block in blocks:
        value *= state[block]
    return value


@lru_cache(maxsize=None)
def local_mass(t: int, ket_count: int, bra_count: int) -> float:
    if ket_count + bra_count == 0:
        return 1.0
    states = occupations(t)
    ket_values = np.array([monomial(state, range(ket_count)) for state in states], dtype=float)
    bra_values = np.array(
        [monomial(state, range(ket_count, ket_count + bra_count)) for state in states],
        dtype=float,
    )
    weights = cp.Variable(len(states), nonneg=True)
    moments = cp.hstack([ket_values @ weights, bra_values @ weights])
    problem = cp.Problem(cp.Maximize(cp.geo_mean(moments)), [cp.sum(weights) == 1])
    value = problem.solve(
        solver="CLARABEL",
        tol_gap_abs=2e-10,
        tol_gap_rel=2e-10,
        tol_feas=2e-10,
        max_iter=500,
    )
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise AssertionError(("local mass", t, ket_count, bra_count, problem.status))
    return float(value)


def integer_partitions(total: int, maximum: int | None = None):
    if total == 0:
        yield ()
        return
    if maximum is None or maximum > total:
        maximum = total
    for first in range(maximum, 0, -1):
        for rest in integer_partitions(total - first, first):
            yield (first,) + rest


def partition_constant(doses: tuple[int, ...]) -> float:
    entries = 2 * len(doses)
    total = 0.0
    for placement in product(range(entries), repeat=4):
        local = 1.0
        maximum_entry = 0
        for node, dose in enumerate(doses):
            ket_count = sum(entry == 2 * node for entry in placement)
            bra_count = sum(entry == 2 * node + 1 for entry in placement)
            maximum_entry = max(maximum_entry, ket_count, bra_count)
            local *= local_mass(dose, ket_count, bra_count)
        if maximum_entry <= 1:
            relative_graph = 1.0
        else:
            # Graph factor N^{(1-k)/2}, divided by the target N^{-1/2}.
            relative_graph = N ** ((2 - maximum_entry) / 2)
        total += relative_graph * local
    return total


def main() -> None:
    local_rows = []
    for t in range(1, DOSE + 1):
        values = []
        for a in range(5):
            for b in range(5 - a):
                if a + b:
                    values.append(f"({a},{b})={local_mass(t,a,b):.8g}")
        local_rows.append(f"t={t}:" + ",".join(values))

    results = []
    for used in range(1, DOSE + 1):
        for doses in integer_partitions(used):
            results.append((partition_constant(doses), doses))
    results.sort(reverse=True)
    worst, worst_partition = results[0]
    print("local block masses:\n" + "\n".join(local_rows))
    print("largest adaptive minimal-chain ledgers:")
    for value, doses in results[:12]:
        print(f"doses={doses},constant={value:.12g},TV_at_N1024={value/32:.12g}")
    if worst_partition != (1, 1, 1, 1, 1, 1) or not np.isclose(worst, 8730.0, atol=2e-4):
        raise AssertionError(("coarse-ledger reference changed", worst_partition, worst))
    if worst <= 32 / 3:
        raise AssertionError(("coarse ledger unexpectedly met threshold", worst))
    print(
        "termwise adaptive ledger rejected: "
        f"worst_partition={worst_partition}, constant={worst:.12g}, "
        f"threshold_overshoot={worst/(32/3):.12g}"
    )


if __name__ == "__main__":
    main()
