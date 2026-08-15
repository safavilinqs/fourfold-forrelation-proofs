#!/usr/bin/env python3
"""Exact rejection of an unweighted l2 temporal-placement ledger."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
import math


ENTRIES = 12  # ket/bra entries of six dose-one nodes


def placement_weight(placement: tuple[int, ...]) -> Fraction:
    value = Fraction(1)
    for node in range(6):
        ket = sum(entry == 2 * node for entry in placement)
        bra = sum(entry == 2 * node + 1 for entry in placement)
        if ket > 1 or bra > 1:
            return Fraction(0)
        if ket == 1 and bra == 1:
            value *= Fraction(1, 2)
    return value


def main() -> None:
    weights = [placement_weight(placement) for placement in product(range(ENTRIES), repeat=4)]
    l1 = sum(weights, Fraction(0))
    l2_squared = sum((weight * weight for weight in weights), Fraction(0))
    nonzero = sum(weight > 0 for weight in weights)
    expected_l1 = Fraction(8730)
    expected_l2_squared = Fraction(14445, 2)
    if l1 != expected_l1 or l2_squared != expected_l2_squared:
        raise AssertionError(("temporal norm reference", l1, l2_squared))
    l2 = math.sqrt(float(l2_squared))
    threshold = 32 / 3
    if l2 <= threshold:
        raise AssertionError(("l2 barrier disappeared", l2, threshold))
    print(
        "temporal square-function barrier confirmed: "
        f"nonzero={nonzero}, l1={l1}, l2_squared={l2_squared}, "
        f"l2={l2:.12g}, threshold_overshoot={l2/threshold:.12g}"
    )


if __name__ == "__main__":
    main()
