#!/usr/bin/env python3
"""Counterexample to the entrywise signed-permutation recording target.

For one passive batch, uniform weights on two dose-two sectors give TV 1/N
on the adjacent split and 1/sqrt(N) on either crossing split.  The entrywise
three-match budget gives only 2/q^3.  The crossing ratio is q^2/2 and rules
out both an entrywise recording constant and a uniform 1/N target.
"""

from __future__ import annotations

import numpy as np


def sylvester(n: int) -> np.ndarray:
    h = np.array([[1.0]])
    while h.shape[0] < n:
        h = np.block([[h, h], [h, -h]])
    return h / np.sqrt(n)


def numeric_advantages(q: int) -> tuple[float, float]:
    n = q * q
    h = sylvester(n)
    tensor = np.einsum("ij,jk,kl->ijkl", h, h, h)
    adjacent = tensor.reshape(n * n, n * n)
    crossing = tensor.transpose(0, 2, 1, 3).reshape(n * n, n * n)
    adjacent_nuclear = np.linalg.norm(adjacent, ord="nuc")
    crossing_nuclear = np.linalg.norm(crossing, ord="nuc")
    if not np.isclose(adjacent_nuclear, n, atol=2e-8):
        raise AssertionError(("adjacent nuclear norm", q, adjacent_nuclear, n))
    if not np.isclose(crossing_nuclear, n ** 1.5, atol=2e-8):
        raise AssertionError(("crossing nuclear norm", q, crossing_nuclear, n ** 1.5))
    # Moment difference is 2T.  Both sectors have uniform total mass 1/2,
    # and TV contributes one half of the full Hermitian trace norm.
    return float(adjacent_nuclear / (n * n)), float(crossing_nuclear / (n * n))


def main() -> None:
    rows = []
    for q in (2, 4):
        adjacent, crossing = numeric_advantages(q)
        if not np.isclose(adjacent, 1 / (q * q), atol=2e-10):
            raise AssertionError(("adjacent advantage", q, adjacent))
        if not np.isclose(crossing, 1 / q, atol=2e-10):
            raise AssertionError(("crossing advantage", q, crossing))
        rows.append(f"q={q}:adjacent={adjacent:.12g},crossing={crossing:.12g}")

    q = 32
    adjacent = 1 / (q * q)
    crossing = 1 / q
    entrywise_record = 2 / (q**3)
    ratio = crossing / entrywise_record
    if not np.isclose(ratio, q * q / 2, atol=1e-15):
        raise AssertionError(("recording gap", ratio, q * q / 2))
    print(
        "minimal-chain recording target falsified: "
        + ", ".join(rows)
        + f", q=32 adjacent={adjacent:.12g}, crossing={crossing:.12g}, "
        f"entrywise_budget={entrywise_record:.12g}, crossing_gap={ratio:.12g}"
    )


if __name__ == "__main__":
    main()
