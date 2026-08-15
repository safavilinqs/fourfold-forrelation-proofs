#!/usr/bin/env python3
"""Exact-arithmetic checks for the constant ledger in the v3.1 derivation.

This script checks arithmetic, not the analytic lemmas in ``main.tex``.
The two classical inputs 25/8 < pi < 22/7 are stated in the manuscript;
their rational consequences are checked below.  The bound e < 11/4 is
also certified below using a rational upper bound on its power series.
"""
from fractions import Fraction
from math import factorial

# A rational certificate for e < 11/4.  For n >= 8, successive terms in
# the exponential series have ratio at most 1/9, so the tail is bounded by
# (1/8!) / (1 - 1/9).
exp_partial = sum((Fraction(1, factorial(n)) for n in range(8)), Fraction())
exp_upper = exp_partial + Fraction(9, 8 * factorial(8))
assert exp_upper == Fraction(876_809, 322_560)
assert exp_upper < Fraction(11, 4)

# Rational consequences of the classical bounds 25/8 < pi < 22/7.
rho = Fraction(200, 201)
beta_lower = 2 * Fraction(7, 22) * rho
assert beta_lower == Fraction(1400, 2211)
assert Fraction(400, 1) / Fraction(25, 8) == 128

# Mean margin: (1400/2211)^3 - 1/4 > 1/260.
margin = beta_lower**3 - Fraction(1, 4)
assert margin == Fraction(167_480_069, 43_234_079_724)
assert margin > Fraction(1, 260)

# Variance/conditioning error at N0 = 2^30.
variance_constant = 1540 * 260**2
assert variance_constant == 104_104_000
n0 = 2**30
conditioning_error = Fraction(variance_constant, n0)
assert conditioning_error == Fraction(1_626_625, 16_777_216)

# Uniform average-case bound for c0 = 2/15.
# e < 11/4 and 8^(1/4) > 5/3 imply q < 11/25.
assert 8 * 3**4 > 5**4
q_upper = Fraction(11, 25)
assert 2 * Fraction(2, 15) * Fraction(11, 4) * Fraction(3, 5) == q_upper
average_gap_upper = 2 * q_upper**4 / (1 - q_upper)
assert average_gap_upper == Fraction(14_641, 109_375)

# Exact-promise conditioning leaves a larger gap.
hard_gap_lower = Fraction(1, 3) - 2 * conditioning_error
assert hard_gap_lower == Fraction(3_508_733, 25_165_824)
assert hard_gap_lower > average_gap_upper
difference = hard_gap_lower - average_gap_upper
assert difference == Fraction(15_314_842_691, 2_752_512_000_000)

# The theorem is nontrivial at N0: (2/15)*N0^(1/8) > 1.
assert 2**8 * n0 > 15**8

print("margin =", margin, "≈", float(margin))
print("rational upper bound on e =", exp_upper, "≈", float(exp_upper))
print("conditioning error at 2^30 =", conditioning_error, "≈", float(conditioning_error))
print("average gap upper bound =", average_gap_upper, "≈", float(average_gap_upper))
print("hard gap lower bound =", hard_gap_lower, "≈", float(hard_gap_lower))
print("gap difference =", difference, "≈", float(difference))
print("c0 = 2/15, alpha = 1/8, N0 = 2^30")
