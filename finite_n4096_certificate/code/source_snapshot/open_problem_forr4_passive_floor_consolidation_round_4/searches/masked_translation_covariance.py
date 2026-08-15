#!/usr/bin/env python3
"""Exact translation covariance behind the masked-contraction reduction.

For one signed-permutation link, translate the left coordinates by ``(a,b)``
and the right coordinates by ``(c,d)``.  Relabelling the hidden permutation
shows that the moment changes only by an explicit Walsh character.  Across
the three-link chain these characters multiply.  On a physical occurrence
entry the full support is a disjoint union of the row and column supports, so
the character factors into a row sign and a column sign.

Consequently simultaneous translation of both diagonal laws changes the
masked occurrence matrix only by row/column permutations and diagonal signs.
The weighted nuclear norm is invariant, and joint concavity implies that
averaging both laws over the full translation group cannot decrease it.
"""

from __future__ import annotations

from collections.abc import Iterable


Support = tuple[int, ...]
Shift = tuple[int, int]
Configuration = tuple[Support, ...]


def walsh(left: int, right: int) -> int:
    return -1 if (left & right).bit_count() % 2 else 1


def translate_support(
    order: int,
    support: Iterable[int],
    shift: Shift,
) -> Support:
    row_shift, column_shift = shift
    return tuple(
        sorted(
            (((coordinate // order) ^ row_shift) * order)
            + ((coordinate % order) ^ column_shift)
            for coordinate in support
        )
    )


def translate_configuration(
    order: int,
    configuration: Configuration,
    shifts: tuple[Shift, ...],
) -> Configuration:
    if len(configuration) != len(shifts):
        raise ValueError("one translation required per block")
    return tuple(
        translate_support(order, support, shift)
        for support, shift in zip(configuration, shifts, strict=True)
    )


def link_translation_sign(
    order: int,
    left_support: Iterable[int],
    right_support: Iterable[int],
    left_shift: Shift,
    right_shift: Shift,
) -> int:
    """Character acquired by one exact signed-permutation link moment."""

    a, b = left_shift
    c, d = right_shift
    if any(value < 0 or value >= order for value in (a, b, c, d)):
        raise ValueError(("translation outside order", order, left_shift, right_shift))
    sign = 1
    for coordinate in left_support:
        row, column = divmod(coordinate, order)
        sign *= walsh(row, c) * walsh(a, c) * walsh(column, d)
    for coordinate in right_support:
        row, column = divmod(coordinate, order)
        sign *= walsh(a, row) * walsh(b, column) * walsh(b, d)
    return sign


def chain_translation_sign(
    order: int,
    supports: Configuration,
    shifts: tuple[Shift, Shift, Shift, Shift],
) -> int:
    if len(supports) != 4:
        raise ValueError("the chain has four support blocks")
    result = 1
    for index in range(3):
        result *= link_translation_sign(
            order,
            supports[index],
            supports[index + 1],
            shifts[index],
            shifts[index + 1],
        )
    return result


def xor_shifts(
    left: tuple[Shift, Shift, Shift, Shift],
    right: tuple[Shift, Shift, Shift, Shift],
) -> tuple[Shift, Shift, Shift, Shift]:
    return tuple(
        (left_row ^ right_row, left_column ^ right_column)
        for (left_row, left_column), (right_row, right_column) in zip(
            left, right, strict=True
        )
    )


def translation_cocycle(
    order: int,
    block_sizes: tuple[int, int, int, int],
    left: tuple[Shift, Shift, Shift, Shift],
    right: tuple[Shift, Shift, Shift, Shift],
) -> int:
    """Scalar projective cocycle on one fixed occurrence-size sector."""

    canonical = tuple(tuple(range(size)) for size in block_sizes)
    combined = xor_shifts(left, right)
    return (
        chain_translation_sign(order, canonical, right)
        * chain_translation_sign(
            order,
            translate_configuration(order, canonical, right),
            left,
        )
        * chain_translation_sign(order, canonical, combined)
    )


def disjoint_union(
    row: Configuration,
    column: Configuration,
) -> Configuration:
    if len(row) != len(column):
        raise ValueError("row and column must have the same block count")
    result = []
    for left, right in zip(row, column, strict=True):
        if set(left).intersection(right):
            raise ValueError("physical occurrence supports must be disjoint")
        result.append(tuple(sorted(left + right)))
    return tuple(result)


def separated_occurrence_signs(
    order: int,
    row: Configuration,
    column: Configuration,
    shifts: tuple[Shift, Shift, Shift, Shift],
) -> tuple[int, int, int]:
    """Return row, column, and full translation characters."""

    full = disjoint_union(row, column)
    row_sign = chain_translation_sign(order, row, shifts)
    column_sign = chain_translation_sign(order, column, shifts)
    full_sign = chain_translation_sign(order, full, shifts)
    if row_sign * column_sign != full_sign:
        raise AssertionError("translation character failed to separate")
    return row_sign, column_sign, full_sign
