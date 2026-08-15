#!/usr/bin/env python3
"""Exact finite occurrence budgets for the signed-permutation hard plant."""

from __future__ import annotations

from fractions import Fraction
from math import comb


def link_budget(q: int, left: int, right: int) -> Fraction:
    return sum(
        (
            Fraction(comb(left, r) * comb(right, r), comb(q, r))
            for r in range(1, min(left, right, q) + 1)
        ),
        Fraction(),
    )


def maximum_budget(q: int, dose: int) -> tuple[Fraction, tuple[int, int, int, int]]:
    best = Fraction()
    allocation = (0, 0, 0, 0)
    total = 2 * dose
    for n1 in range(total + 1):
        for n2 in range(total - n1 + 1):
            for n3 in range(total - n1 - n2 + 1):
                for n4 in range(total - n1 - n2 - n3 + 1):
                    value = (
                        link_budget(q, n1, n2)
                        * link_budget(q, n2, n3)
                        * link_budget(q, n3, n4)
                    )
                    if value > best:
                        best = value
                        allocation = (n1, n2, n3, n4)
    return best, allocation


def main() -> None:
    print("q N dose allocation signed_budget")
    for q in (16, 32, 64):
        for dose in range(4, 9):
            budget, allocation = maximum_budget(q, dose)
            print(q, q * q, dose, allocation, f"{float(2 * budget):.12g}")

    target, allocation = maximum_budget(32, 6)
    signed = 2 * target
    slack = Fraction(1, 3) / signed
    expected = Fraction(17_497_415, 442_336_768)
    if target != expected:
        raise AssertionError(("N=1024 exact budget changed", target, expected))
    if signed >= Fraction(1, 3):
        raise AssertionError(("dose-six budget misses threshold", signed))
    print(
        "N=1024 dose=6 certificate: "
        f"allocation={allocation}, one_sided={target}, "
        f"signed={float(signed):.12g}, admissible_extra_factor={float(slack):.12g}"
    )


if __name__ == "__main__":
    main()
