#!/usr/bin/env python3
"""Exact Fourier moments of one signed-permutation plant link.

For a uniformly random signed permutation ``P`` of order ``q``, let

``x = (K P).reshape(-1)`` and ``y = (P K).reshape(-1)``,

where ``K`` is the unnormalised Sylvester matrix. This module evaluates

``E[prod(x[i] for i in left) prod(y[j] for j in right)]``

exactly, without enumerating the ``q! 2^q`` signed permutations. Only the
rows and columns touched by the two small Fourier supports enter the
calculation, so the dose-six use case remains inexpensive at ``q=32``.

The sign average first forces the odd column degrees of ``left`` to match
the odd row degrees of ``right`` under the hidden permutation. The remaining
permutation average is a permanent. Its odd block is small, and its even
block differs from the all-ones matrix only in touched rows and columns. A
special-row/special-column expansion evaluates that permanent exactly with
integer arithmetic.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from fractions import Fraction
from functools import lru_cache
from math import factorial


Support = tuple[int, ...]


def _walsh(row: int, column: int) -> int:
    """Return one entry of the unnormalised Sylvester matrix."""

    return -1 if (row & column).bit_count() % 2 else 1


def _permanent(matrix: tuple[tuple[int, ...], ...]) -> int:
    """Permanent of a small integer matrix by subset dynamic programming."""

    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("permanent requires a square matrix")
    dynamic = [0] * (1 << size)
    dynamic[0] = 1
    for row_index, row in enumerate(matrix):
        updated = [0] * (1 << size)
        for mask, value in enumerate(dynamic):
            if not value or mask.bit_count() != row_index:
                continue
            for column, entry in enumerate(row):
                if not (mask >> column) & 1:
                    updated[mask | (1 << column)] += value * entry
        dynamic = updated
    return dynamic[-1]


def _injection_sum(
    sources: tuple[int, ...],
    targets: tuple[int, ...],
    weight: Callable[[int, int], int],
) -> int:
    """Sum products of ``weight`` over all injections source -> target."""

    if len(sources) > len(targets):
        return 0
    dynamic = [0] * (1 << len(sources))
    dynamic[0] = 1
    for target in targets:
        updated = dynamic.copy()
        for mask, value in enumerate(dynamic):
            if not value:
                continue
            for index, source in enumerate(sources):
                if not (mask >> index) & 1:
                    updated[mask | (1 << index)] += value * weight(source, target)
        dynamic = updated
    return dynamic[-1]


def _indices(mask: int, values: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(value for index, value in enumerate(values) if (mask >> index) & 1)


def _structured_permanent(
    domain: tuple[int, ...],
    codomain: tuple[int, ...],
    left_xors: dict[int, int],
    right_xors: dict[int, int],
) -> int:
    """Permanent of the even-parity block.

    The entry at domain column ``b`` and codomain row ``c`` is
    ``K[left_xors[b], c] K[b, right_xors[c]]``. It is one unless ``b`` or
    ``c`` is touched nontrivially. The expansion partitions permutations by
    which special domain elements map to which special codomain elements.
    """

    if len(domain) != len(codomain):
        raise ValueError("permanent blocks must be square")
    special_domain = tuple(value for value in domain if left_xors.get(value, 0))
    special_codomain = tuple(value for value in codomain if right_xors.get(value, 0))
    special_domain_set = set(special_domain)
    special_codomain_set = set(special_codomain)
    ordinary_domain = tuple(
        value for value in domain if value not in special_domain_set
    )
    ordinary_codomain = tuple(
        value for value in codomain if value not in special_codomain_set
    )

    def left_weight(source: int, target: int) -> int:
        return _walsh(left_xors[source], target)

    def right_weight(source: int, target: int) -> int:
        # The source is a special codomain row and the target is its
        # preimage in the permutation domain.
        return _walsh(target, right_xors[source])

    left_injections: dict[int, int] = {}
    for selected_mask in range(1 << len(special_domain)):
        remaining = tuple(
            value
            for index, value in enumerate(special_domain)
            if not (selected_mask >> index) & 1
        )
        left_injections[selected_mask] = _injection_sum(
            remaining, ordinary_codomain, left_weight
        )

    right_injections: dict[int, int] = {}
    for selected_mask in range(1 << len(special_codomain)):
        remaining = tuple(
            value
            for index, value in enumerate(special_codomain)
            if not (selected_mask >> index) & 1
        )
        right_injections[selected_mask] = _injection_sum(
            remaining, ordinary_domain, right_weight
        )

    total = 0
    domain_size = len(special_domain)
    codomain_size = len(special_codomain)
    block_size = len(domain)
    for domain_mask in range(1 << domain_size):
        shared = domain_mask.bit_count()
        for codomain_mask in range(1 << codomain_size):
            if codomain_mask.bit_count() != shared:
                continue
            selected_domain = _indices(domain_mask, special_domain)
            selected_codomain = _indices(codomain_mask, special_codomain)
            joint = _permanent(
                tuple(
                    tuple(
                        left_weight(source, target) * right_weight(target, source)
                        for target in selected_codomain
                    )
                    for source in selected_domain
                )
            )
            remaining = block_size - domain_size - codomain_size + shared
            if remaining < 0:
                # Too few ordinary elements remain to realise this choice
                # of special-to-special matches.
                continue
            total += (
                joint
                * left_injections[domain_mask]
                * right_injections[codomain_mask]
                * factorial(remaining)
            )
    return total


def _validate_support(order: int, support: Iterable[int]) -> Support:
    result = tuple(sorted(support))
    if len(result) != len(set(result)):
        raise ValueError("Fourier supports cannot contain duplicate coordinates")
    if any(coordinate < 0 or coordinate >= order * order for coordinate in result):
        raise ValueError(("coordinate outside q^2 block", order, result))
    return result


@lru_cache(maxsize=None)
def link_moment(
    order: int,
    left_support: Support,
    right_support: Support,
) -> Fraction:
    """Return the exact signed-permutation link moment.

    Supports must be sorted tuples of distinct flattened coordinates. Use
    :func:`moment` for arbitrary iterables.
    """

    left = _validate_support(order, left_support)
    right = _validate_support(order, right_support)

    left_counts: dict[int, int] = {}
    left_xors: dict[int, int] = {}
    for coordinate in left:
        row, column = divmod(coordinate, order)
        left_counts[column] = left_counts.get(column, 0) + 1
        left_xors[column] = left_xors.get(column, 0) ^ row

    right_counts: dict[int, int] = {}
    right_xors: dict[int, int] = {}
    for coordinate in right:
        row, column = divmod(coordinate, order)
        right_counts[row] = right_counts.get(row, 0) + 1
        right_xors[row] = right_xors.get(row, 0) ^ column

    odd_domain = tuple(value for value in range(order) if left_counts.get(value, 0) % 2)
    odd_codomain = tuple(
        value for value in range(order) if right_counts.get(value, 0) % 2
    )
    if len(odd_domain) != len(odd_codomain):
        return Fraction(0)

    odd_matrix = tuple(
        tuple(
            _walsh(left_xors.get(domain, 0), codomain)
            * _walsh(domain, right_xors.get(codomain, 0))
            for codomain in odd_codomain
        )
        for domain in odd_domain
    )
    odd_permanent = _permanent(odd_matrix)
    if not odd_permanent:
        return Fraction(0)

    odd_domain_set = set(odd_domain)
    odd_codomain_set = set(odd_codomain)
    even_domain = tuple(value for value in range(order) if value not in odd_domain_set)
    even_codomain = tuple(
        value for value in range(order) if value not in odd_codomain_set
    )
    even_permanent = _structured_permanent(
        even_domain, even_codomain, left_xors, right_xors
    )
    return Fraction(odd_permanent * even_permanent, factorial(order))


def moment(
    order: int,
    left_support: Iterable[int],
    right_support: Iterable[int],
) -> Fraction:
    """Validated iterable-friendly wrapper around :func:`link_moment`."""

    left = _validate_support(order, left_support)
    right = _validate_support(order, right_support)
    return link_moment(order, left, right)


def chain_moment(order: int, supports: tuple[Iterable[int], ...]) -> Fraction:
    """Exact four-block plus-plant moment for three independent links."""

    if len(supports) != 4:
        raise ValueError("the chain plant has four blocks")
    canonical = tuple(_validate_support(order, support) for support in supports)
    result = Fraction(1)
    for left, right in zip(canonical, canonical[1:]):
        result *= link_moment(order, left, right)
        if not result:
            break
    return result
