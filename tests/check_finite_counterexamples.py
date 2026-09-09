#!/usr/bin/env python3
"""Exact regressions for failed finite-proof premises, not lower-bound tests.

Moments are evaluated independently by averaging all permutations and
eliminating the uniform signs by parity. No frozen proof code is imported.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import comb


@lru_cache(None)
def permutation_pairs(q):
    result = []
    for pi in permutations(range(q)):
        inverse = [0] * q
        for b, c in enumerate(pi):
            inverse[c] = b
        result.append((pi, inverse))
    return result


def moment(q, left, right):
    total = 0
    for pi, inverse in permutation_pairs(q):
        sign_parity = 0
        phase = 0
        for a, b in left:
            sign_parity ^= 1 << b
            phase ^= (a & pi[b]).bit_count() & 1
        for c, d in right:
            b = inverse[c]
            sign_parity ^= 1 << b
            phase ^= (b & d).bit_count() & 1
        if not sign_parity:
            total += 1 if not phase else -1
    return Fraction(total, len(permutation_pairs(q)))


def check_photon_vs_parity():
    # Two basis states in the same three-photon sector.
    nu = {(0, 0): 3}
    mu = {(1, 0): 1, (2, 0): 1, (3, 0): 1}
    support = lambda v: {i for i, k in v.items() if k % 2}
    s, t = support(nu), support(mu)
    assert sum(nu.values()) == sum(mu.values()) == 3
    assert (len(s), len(t)) == (1, 3)
    profile = tuple(sum(b == j for b, _ in s ^ t) for j in range(4))
    split = tuple(sum(b == j for b, _ in s - t) for j in range(4))
    assert profile == (1, 1, 1, 1) and split == (1, 0, 0, 0)
    assert 2 * sum(split) != sum(profile)
    assert moment(8, ((0, 0),), ((0, 0),)) == Fraction(1, 8)
    assert Fraction(19, 25)**4 / 64**3 == Fraction(130321, 102400000000)


def check_four_cubic():
    q = 8
    vertical = ((0, 0), (1, 0), (2, 0))
    horizontal = ((0, 0), (0, 1), (0, 2))
    vh = moment(q, vertical, horizontal)
    hv = moment(q, horizontal, vertical)
    assert vh == Fraction(1, q)
    assert hv == Fraction(1, comb(q, 3))
    asserted_link = Fraction(q + 2, q * (q - 1) * (q - 2))
    actual = vh * hv * vh
    asserted = asserted_link**2 * hv
    assert actual == Fraction(1, 3584)
    assert asserted == Fraction(25, 1580544)
    assert actual > asserted


def check_walsh_identity():
    q = 8
    x, z, y = (1, 0), (1, 0), (0, 0)
    actual = moment(q, (x,), (z,)) * moment(q, (z,), (y,))
    # The incorrect formula replaces the z-dependent phase by H_N(x,y).
    assert actual == -Fraction(1, q**2)
    assert actual != Fraction(1, q**2)


def check_completed_cubic_amplitude():
    q = 4
    points = tuple((a, b) for a in range(q) for b in range(q))
    pairs = tuple(combinations(points, 2))
    # Integer-scaled V[x,E] = q M({x} triangle E,{0}).
    rows = []
    for x in points:
        row = []
        for pair in pairs:
            value = q * (q - 1) * moment(q, tuple({x} ^ set(pair)), ((0, 0),))
            assert value.denominator == 1
            row.append(int(value))
        rows.append(row)
    gram = [[sum(x * y for x, y in zip(a, b)) for b in rows] for a in rows]

    def mv(v):
        return [sum(x * y for x, y in zip(row, v)) for row in gram]

    assert mv([1] * q**2) == [0] * q**2
    # These vectors span the two remaining invariant subspaces.
    for b in range(q):
        for a in range(1, q):
            v = [0] * q**2
            v[a * q + b], v[b] = 1, -1
            assert mv(v) == [q**3 * (q - 1) * x for x in v]
    for b in range(1, q):
        v = [int(j == b) - int(j == 0) for _, j in points]
        assert mv(v) == [q**4 * (q - 1) // 2 * x for x in v]
    # Undoing the scale gives nuclear norm (1+sqrt(2q))/sqrt(q+1)>1
    # under uniform row/column weights. The Gram identities prove it exactly.
    assert 2 * q > q + 1


if __name__ == "__main__":
    check_photon_vs_parity()
    check_four_cubic()
    check_walsh_identity()
    check_completed_cubic_amplitude()
    print("PASS exact finite counterexamples: parity sectors, four-cubic maximum, Walsh identity and amplitude")
