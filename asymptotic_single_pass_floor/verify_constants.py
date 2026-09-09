#!/usr/bin/env python3
"""Exact-arithmetic checks for the parallel-query derivation.

This script checks arithmetic, not the analytic lemmas in ``main.tex``.
The two classical inputs 25/8 < pi < 22/7 are stated in the manuscript;
their rational consequences are checked below.  The bound e < 3 is
also certified below using a rational upper bound on its power series.
"""
from fractions import Fraction

# For n >= 2, successive terms in the exponential series have ratio at
# most 1/3, so e <= 1 + 1 + (1/2!) / (1 - 1/3) = 11/4 < 3.
exp_upper = 2 + Fraction(1, 2) / (1 - Fraction(1, 3))
assert exp_upper == Fraction(11, 4)
assert exp_upper < 3

# Rational consequences of the classical bounds 25/8 < pi < 22/7.
rho = Fraction(200, 201)
beta_lower = 2 * Fraction(7, 22) * rho
assert beta_lower == Fraction(1400, 2211)
assert Fraction(400, 1) / Fraction(25, 8) == 128

# Mean margin: (1400/2211)^3 - 1/4 > 1/260.
margin = beta_lower**3 - Fraction(1, 4)
assert margin == Fraction(167_480_069, 43_234_079_724)
assert margin > Fraction(1, 260)

# Variance and promise-failure probability at N0 = 2^30.
variance_prefactor = 4 + 12 * 128
assert variance_prefactor == 1540
variance_constant = variance_prefactor * 260**2
assert variance_constant == 104_104_000
n0 = 2**30
promise_failure = Fraction(variance_constant, n0)
assert promise_failure == Fraction(1_626_625, 16_777_216)
assert promise_failure < Fraction(1, 10)

# Uniform average-case bound for c0 = 2/15.
# e < 3 and 2^(1/4) < 5/4 imply q = 2^(1/4)e D N^(-1/8) < 1/2.
assert 2 < Fraction(5, 4)**4
q_upper = Fraction(2, 15) * 3 * Fraction(5, 4)
assert q_upper == Fraction(1, 2)
parallel_gap_upper = q_upper**4 / (1 - q_upper**2)
assert parallel_gap_upper == Fraction(1, 12)

# The sufficient adaptive Fourier bound is unhalved; it is not proved by
# this arithmetic check. Its even-level series would still close the gap.
adaptive_gap_upper = 2 * parallel_gap_upper
assert adaptive_gap_upper == Fraction(1, 6)

# Direct use of the two promise masses leaves a gap greater than 1/5.
hard_gap_lower = Fraction(1, 3) - Fraction(4, 3) * promise_failure
hard_gap_floor = Fraction(1, 3) - Fraction(4, 3) * Fraction(1, 10)
assert hard_gap_floor == Fraction(1, 5)
assert parallel_gap_upper < adaptive_gap_upper < hard_gap_floor < hard_gap_lower

# The theorem is nontrivial at N0: (2/15)*N0^(1/8) > 1.
assert 2**8 * n0 > 15**8

print("margin =", margin, "≈", float(margin))
print("rational upper bound on e =", exp_upper, "≈", float(exp_upper))
print("promise-failure bound at 2^30 =", promise_failure, "≈", float(promise_failure))
print("parallel gap upper bound =", parallel_gap_upper)
print("conditional adaptive gap upper bound =", adaptive_gap_upper)
print("hard gap lower bound =", hard_gap_lower, "≈", float(hard_gap_lower))
print("closure: 1/12 < 1/6 < 1/5 < hard gap")
print("c0 = 2/15, alpha = 1/8, N0 = 2^30")
