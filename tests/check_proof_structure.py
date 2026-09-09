#!/usr/bin/env python3
"""Exact checks of the record identity and a counterexample to the norm inference."""
from fractions import Fraction
from itertools import combinations, product


def parity(word):
    return frozenset(x for x in set(word) if x and word.count(x) % 2)


def admissible(word, positions):
    labels = tuple(word[i] for i in positions)
    if 0 in labels or len(set(labels)) != len(labels):
        return None
    if any(word[i] in word[i + 1:] for i in positions):
        return None
    marked = frozenset(labels)
    residual = parity(word) ^ marked
    if residual & marked:
        return None
    return marked, residual


def check_records():
    # Include repeated labels and vacuum on both sides.
    modes, dose = 4, 3
    words = list(product(range(modes + 1), repeat=dose))
    checked = 0
    for ket, bra in product(words, repeat=2):
        difference = parity(ket) ^ parity(bra)
        expected = frozenset(
            (0 if label in parity(ket) else 1,
             max(i for i, value in enumerate(
                 ket if label in parity(ket) else bra) if value == label))
            for label in difference
        )
        accepted = []
        for record in combinations(range(2 * dose), len(difference)):
            left = tuple(i for i in record if i < dose)
            right = tuple(i - dose for i in record if i >= dose)
            a, b = admissible(ket, left), admissible(bra, right)
            if a is None or b is None or a[0] & b[0] or a[1] != b[1]:
                continue
            accepted.append(frozenset(
                [(0, i) for i in left] + [(1, i) for i in right]
            ))
        assert accepted == [expected], (ket, bra, accepted, expected)
        checked += 1
    print(f"PASS exact last-occurrence partition: {checked} word pairs")


def check_norm_counterexample():
    # K=J has gamma_2=1; D=J/n has normalized scalar Gram features.
    # The weighted kernel J/n is PSD rank one with trace norm one.
    for n in (2, 3, 16):
        feature_mass = sum([Fraction(1, n)] * n)
        weighted_trace_norm = sum([Fraction(1, n)] * n)
        entry_sum = sum([Fraction(1, n)] * (n * n))
        assert feature_mass == weighted_trace_norm == 1
        assert entry_sum == n > weighted_trace_norm
    print("PASS counterexample: feature normalization does not bound entry sum")


if __name__ == "__main__":
    check_records()
    check_norm_counterexample()
