#!/usr/bin/env python3
"""Direct small-order witness for the leading adjacent cubic--quintic split.

For profile ``(1,1,3,5)`` and split ``(0,1,2,2)``, rows are indexed by

    (b, C, D)

and columns by

    (a, e, T).

Here ``a,b,e`` are singleton cells, ``C,D`` are selected endpoint pairs,
and ``T`` is the remaining quintic triple.  The full supports are
``C union {e}`` and ``D union T``.  The matrix entry is

    M_11(a,b) M_13(b,C union {e}) M_35(C union {e},D union T).

The complementary split is its transpose.  This script computes the link
moments exactly by averaging every signed permutation at q=4, restricts the
two pairs and the triple to translation orbits, and evaluates the normalized
nuclear norm.  It is a physical lower-witness diagnostic, not a q=32 claim.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, sqrt

import numpy as np
from scipy import sparse

from opposite_endpoint_orbit_scan import (
    pair_orbit,
    translated_sequence,
)


@dataclass(frozen=True)
class AdjacentOrbitWitness:
    order: int
    cubic_pair_difference: int
    quintic_pair_difference: int
    triple_shape: tuple[int, int, int]
    rows: int
    columns: int
    nonzero: int
    rank: int
    coefficient: float
    normalized_operator: float


@dataclass(frozen=True)
class ReducedAdjacentOrbitWitness:
    order: int
    cubic_pair_difference: int
    quintic_pair_difference: int
    triple_shape: tuple[int, int, int]
    block_rows: int
    block_columns: int
    nonzero: int
    rank: int
    coefficient: float
    normalized_operator: float


@dataclass(frozen=True)
class LinkMoments:
    supports_three: tuple[tuple[int, ...], ...]
    supports_five: tuple[tuple[int, ...], ...]
    moment_11: np.ndarray
    moment_13: np.ndarray
    moment_35: np.ndarray


@dataclass(frozen=True)
class AdjacentJointSlice:
    order: int
    fixed_singleton: int
    fixed_triple: tuple[int, int, int]
    cubic_extensions: int
    quintic_extensions: int
    descriptor_counts: tuple[int, int, int]
    maximum_per_descriptor: tuple[float, float, float]
    squared_moment_sum: float
    chained_coefficient: float


@dataclass(frozen=True)
class AdjacentRecordThreeSlice:
    order: int
    fixed_singleton: int
    fixed_triple: tuple[int, int, int]
    cubic_extensions: int
    quintic_extensions: int
    squared_moment_sum: float
    chained_coefficient: float


@dataclass(frozen=True)
class HorizontalAdjacentSliceCertificate:
    order: int
    record_one_extensions: tuple[int, int]
    record_one_squared_sum: float
    record_three_squared_sum: float
    record_one_coefficient: float
    record_three_coefficient: float
    combined_coefficient: float


def unnormalized_sylvester(order: int) -> np.ndarray:
    result = np.asarray([[1]], dtype=np.int8)
    while len(result) < order:
        result = np.block([[result, result], [result, -result]])
    if len(result) != order:
        raise ValueError(("power-of-two order required", order))
    return result


@lru_cache(maxsize=None)
def exact_link_moments(order: int) -> LinkMoments:
    """Average all signed permutations and return the three required links."""

    if order != 4:
        raise ValueError("direct signed-permutation enumeration is calibrated at q=4")
    dimension = order * order
    hadamard = unnormalized_sylvester(order)
    left_values = []
    right_values = []
    for permutation in permutations(range(order)):
        for signs in product((-1, 1), repeat=order):
            signed_permutation = np.zeros((order, order), dtype=np.int8)
            for column, row in enumerate(permutation):
                signed_permutation[row, column] = signs[column]
            left_values.append((hadamard @ signed_permutation).reshape(-1))
            right_values.append((signed_permutation @ hadamard).reshape(-1))
    left = np.asarray(left_values, dtype=np.int8)
    right = np.asarray(right_values, dtype=np.int8)
    supports_three = tuple(combinations(range(dimension), 3))
    supports_five = tuple(combinations(range(dimension), 5))

    def features(
        values: np.ndarray,
        supports: tuple[tuple[int, ...], ...],
    ) -> np.ndarray:
        return np.asarray(
            [np.prod(values[:, support], axis=1) for support in supports],
            dtype=np.int8,
        ).T

    left_one = left
    right_one = right
    right_three = features(right, supports_three)
    left_three = features(left, supports_three)
    right_five = features(right, supports_five)
    normalization = len(left)
    return LinkMoments(
        supports_three=supports_three,
        supports_five=supports_five,
        moment_11=left_one.T.astype(float) @ right_one / normalization,
        moment_13=left_one.T.astype(float) @ right_three / normalization,
        moment_35=left_three.T.astype(float) @ right_five / normalization,
    )


def character(left: int, right: int) -> int:
    return -1 if int(left & right).bit_count() % 2 else 1


def xor_values(values: list[int]) -> int:
    result = 0
    for value in values:
        result ^= value
    return result


def character_sum_excluding(
    order: int,
    label: int,
    excluded: tuple[int, ...],
) -> int:
    """Sum a Walsh character over F_q with the named points removed."""

    if label == 0:
        return order - len(excluded)
    return -sum(character(label, value) for value in excluded)


def record_one_link_moment(
    order: int,
    left_support: tuple[int, ...],
    right_support: tuple[int, ...],
) -> float:
    """Exact signed-permutation moment when both parity records have size one.

    After matching the unique odd left column to the unique odd right row,
    a degree-three left support has at most one even column group and a
    degree-five right support has at most two even row groups.  The formulas
    below sum the remaining partial bijection exactly.
    """

    left_columns: dict[int, list[int]] = {}
    for coordinate in left_support:
        row, column = divmod(coordinate, order)
        left_columns.setdefault(column, []).append(row)
    right_rows: dict[int, list[int]] = {}
    for coordinate in right_support:
        row, column = divmod(coordinate, order)
        right_rows.setdefault(row, []).append(column)
    odd_left = [
        column
        for column, rows in left_columns.items()
        if len(rows) % 2
    ]
    odd_right = [
        row
        for row, columns in right_rows.items()
        if len(columns) % 2
    ]
    if len(odd_left) != 1 or len(odd_right) != 1:
        return 0.0
    odd_column = odd_left[0]
    odd_row = odd_right[0]
    odd_left_xor = xor_values(left_columns[odd_column])
    odd_right_xor = xor_values(right_rows[odd_row])
    even_left = [
        (column, xor_values(rows))
        for column, rows in left_columns.items()
        if column != odd_column and xor_values(rows) != 0
    ]
    even_right = [
        (row, xor_values(columns))
        for row, columns in right_rows.items()
        if row != odd_row and xor_values(columns) != 0
    ]
    if len(even_left) > 1 or len(even_right) > 2:
        raise ValueError(
            ("record-one formula degree exceeded", left_support, right_support)
        )

    n = order - 1
    phase = (
        character(odd_left_xor, odd_row)
        * character(odd_column, odd_right_xor)
    )
    if not even_left:
        if not even_right:
            conditional = 1.0
        elif len(even_right) == 1:
            _, label = even_right[0]
            conditional = (
                character_sum_excluding(
                    order,
                    label,
                    (odd_column,),
                )
                / n
            )
        else:
            (_, first_label), (_, second_label) = even_right
            first_sum = character_sum_excluding(
                order,
                first_label,
                (odd_column,),
            )
            second_sum = character_sum_excluding(
                order,
                second_label,
                (odd_column,),
            )
            product_sum = character_sum_excluding(
                order,
                first_label ^ second_label,
                (odd_column,),
            )
            conditional = (
                first_sum * second_sum - product_sum
            ) / (n * (n - 1))
    else:
        even_column, left_label = even_left[0]
        if not even_right:
            conditional = (
                character_sum_excluding(
                    order,
                    left_label,
                    (odd_row,),
                )
                / n
            )
        elif len(even_right) == 1:
            even_row, right_label = even_right[0]
            left_at_row = character(left_label, even_row)
            right_at_column = character(right_label, even_column)
            left_remainder = character_sum_excluding(
                order,
                left_label,
                (odd_row, even_row),
            )
            right_remainder = character_sum_excluding(
                order,
                right_label,
                (odd_column, even_column),
            )
            conditional = (
                left_at_row * right_at_column / n
                + left_remainder
                * right_remainder
                / (n * (n - 1))
            )
        else:
            (
                (first_row, first_label),
                (second_row, second_label),
            ) = even_right
            first_remainder = character_sum_excluding(
                order,
                first_label,
                (odd_column, even_column),
            )
            second_remainder = character_sum_excluding(
                order,
                second_label,
                (odd_column, even_column),
            )
            product_remainder = character_sum_excluding(
                order,
                first_label ^ second_label,
                (odd_column, even_column),
            )
            right_pair_average = (
                first_remainder * second_remainder
                - product_remainder
            ) / ((n - 1) * (n - 2))
            left_remainder = character_sum_excluding(
                order,
                left_label,
                (odd_row, first_row, second_row),
            )
            conditional = (
                character(left_label, first_row)
                * character(first_label, even_column)
                * second_remainder
                / (n * (n - 1))
                + character(left_label, second_row)
                * character(second_label, even_column)
                * first_remainder
                / (n * (n - 1))
                + left_remainder * right_pair_average / n
            )
    return phase * conditional / order


def record_three_link_moment(
    order: int,
    left_support: tuple[int, ...],
    right_support: tuple[int, ...],
) -> float:
    """Exact M_35 entry when both signed-permutation records have size three."""

    left_columns: dict[int, list[int]] = {}
    for coordinate in left_support:
        row, column = divmod(coordinate, order)
        left_columns.setdefault(column, []).append(row)
    right_rows: dict[int, list[int]] = {}
    for coordinate in right_support:
        row, column = divmod(coordinate, order)
        right_rows.setdefault(row, []).append(column)
    odd_columns = tuple(
        sorted(
            column
            for column, rows in left_columns.items()
            if len(rows) % 2
        )
    )
    odd_rows = tuple(
        sorted(
            row
            for row, columns in right_rows.items()
            if len(columns) % 2
        )
    )
    if len(odd_columns) != 3 or len(odd_rows) != 3:
        return 0.0
    left_xors = {
        column: xor_values(left_columns[column])
        for column in odd_columns
    }
    right_xors = {
        row: xor_values(right_rows[row])
        for row in odd_rows
    }
    even_rows = [
        (row, xor_values(columns))
        for row, columns in right_rows.items()
        if row not in odd_rows and xor_values(columns) != 0
    ]
    if len(even_rows) > 1:
        raise ValueError(("record-three quintic degree exceeded", right_support))
    total = 0.0
    for mapped_rows in permutations(odd_rows):
        phase = 1
        for column, row in zip(odd_columns, mapped_rows, strict=True):
            phase *= character(left_xors[column], row)
            phase *= character(column, right_xors[row])
        if even_rows:
            _, label = even_rows[0]
            phase *= (
                -sum(character(label, column) for column in odd_columns)
                / (order - 3)
            )
        total += phase
    return total / (comb(order, 3) * 6)


def parity_descriptor(
    order: int,
    support: tuple[int, ...],
    *,
    axis: int,
) -> tuple[int, int, tuple[tuple[int, int], ...]] | None:
    """Return the odd label/xor and nontrivial even label/xor groups."""

    groups: dict[int, list[int]] = {}
    for coordinate in support:
        row, column = divmod(coordinate, order)
        labels = (row, column)
        groups.setdefault(labels[axis], []).append(labels[1 - axis])
    odd = [
        label for label, values in groups.items() if len(values) % 2
    ]
    if len(odd) != 1:
        return None
    odd_label = odd[0]
    even = tuple(
        sorted(
            (label, xor_values(values))
            for label, values in groups.items()
            if label != odd_label and xor_values(values) != 0
        )
    )
    return odd_label, xor_values(groups[odd_label]), even


def l_shapes_containing(
    order: int,
    fixed_cell: int,
) -> tuple[tuple[int, int, int], ...]:
    """Enumerate all three-edge L-shapes containing one fixed cell."""

    fixed_row, fixed_column = divmod(fixed_cell, order)
    result = set()
    for other_row in range(order):
        if other_row == fixed_row:
            continue
        for other_column in range(order):
            if other_column == fixed_column:
                continue
            # The fixed cell is the corner, horizontal arm, or vertical arm.
            result.add(
                tuple(
                    sorted(
                        (
                            fixed_cell,
                            fixed_row * order + other_column,
                            other_row * order + fixed_column,
                        )
                    )
                )
            )
            result.add(
                tuple(
                    sorted(
                        (
                            fixed_cell,
                            fixed_row * order + other_column,
                            other_row * order + other_column,
                        )
                    )
                )
            )
            result.add(
                tuple(
                    sorted(
                        (
                            fixed_cell,
                            other_row * order + fixed_column,
                            other_row * order + other_column,
                        )
                    )
                )
            )
    expected = 3 * (order - 1) ** 2
    if len(result) != expected:
        raise AssertionError(("L-shape fixed-cell incidence", len(result), expected))
    return tuple(sorted(result))


def endpoint_quintics_containing(
    order: int,
    fixed_triple: tuple[int, int, int],
) -> tuple[tuple[int, ...], ...]:
    """Enumerate quintics containing a triple and having row record one."""

    dimension = order * order
    fixed = set(fixed_triple)
    available = tuple(
        coordinate for coordinate in range(dimension) if coordinate not in fixed
    )
    result = []
    for extra in combinations(available, 2):
        support = tuple(sorted(fixed_triple + extra))
        if parity_descriptor(order, support, axis=0) is not None:
            result.append(support)
    return tuple(result)


def supports_containing_with_records(
    order: int,
    fixed: tuple[int, ...],
    degree: int,
    row_record: int,
    column_record: int | None = None,
) -> tuple[tuple[int, ...], ...]:
    """Enumerate small-codimension extensions with named parity records."""

    dimension = order * order
    fixed_set = set(fixed)
    available = tuple(
        coordinate for coordinate in range(dimension) if coordinate not in fixed_set
    )
    result = []
    for extra in combinations(available, degree - len(fixed)):
        support = tuple(sorted(fixed + extra))
        row = parity_record_size(order, support, axis=0)
        column = parity_record_size(order, support, axis=1)
        if row == row_record and (
            column_record is None or column == column_record
        ):
            result.append(support)
    return tuple(result)


def parity_record_size(
    order: int,
    support: tuple[int, ...],
    *,
    axis: int,
) -> int:
    counts: dict[int, int] = {}
    for coordinate in support:
        labels = divmod(coordinate, order)
        counts[labels[axis]] = counts.get(labels[axis], 0) + 1
    return sum(count % 2 for count in counts.values())


def joint_record_one_slice_sum(
    order: int,
    fixed_singleton: int,
    fixed_triple: tuple[int, int, int],
    *,
    chunk_size: int = 512,
) -> AdjacentJointSlice:
    """Compute the exact M_35 squared slice for one fixed singleton/triple."""

    cubic_supports = l_shapes_containing(order, fixed_singleton)
    quintic_supports = endpoint_quintics_containing(order, fixed_triple)
    left = []
    for support in cubic_supports:
        descriptor = parity_descriptor(order, support, axis=1)
        if descriptor is None or len(descriptor[2]) != 1:
            raise AssertionError(("cubic L descriptor", support, descriptor))
        odd_column, _, ((even_column, left_label),) = descriptor
        left.append((odd_column, even_column, left_label))
    left_array = np.asarray(left, dtype=np.int16)
    odd_column = left_array[:, 0, None]
    even_column = left_array[:, 1, None]
    left_label = left_array[:, 2, None]
    characters = np.asarray(
        [
            [character(label, value) for value in range(order)]
            for label in range(order)
        ],
        dtype=np.int8,
    )
    n = order - 1
    squared_sum = 0.0
    descriptor_counts = [0, 0, 0]
    maximum_per_descriptor = [0.0, 0.0, 0.0]

    for start in range(0, len(quintic_supports), chunk_size):
        selected = quintic_supports[start : start + chunk_size]
        grouped: dict[int, list[tuple[int, ...]]] = {0: [], 1: [], 2: []}
        for support in selected:
            descriptor = parity_descriptor(order, support, axis=0)
            if descriptor is None:
                raise AssertionError(("quintic record descriptor", support))
            odd_row, _, even = descriptor
            if len(even) > 2:
                raise AssertionError(("too many quintic even rows", support))
            flattened = [odd_row]
            for row, label in even:
                flattened.extend((row, label))
            grouped[len(even)].append(tuple(flattened))

        if grouped[0]:
            data = np.asarray(grouped[0], dtype=np.int16)
            odd_row = data[:, 0][None, :]
            conditional = -characters[left_label, odd_row] / n
            by_descriptor = np.square(conditional / order).sum(axis=0)
            squared_sum += float(by_descriptor.sum())
            descriptor_counts[0] += len(grouped[0])
            maximum_per_descriptor[0] = max(
                maximum_per_descriptor[0],
                float(by_descriptor.max()),
            )

        if grouped[1]:
            data = np.asarray(grouped[1], dtype=np.int16)
            odd_row = data[:, 0][None, :]
            even_row = data[:, 1][None, :]
            right_label = data[:, 2][None, :]
            left_at_row = characters[left_label, even_row]
            right_at_column = characters[right_label, even_column]
            left_remainder = -(
                characters[left_label, odd_row]
                + characters[left_label, even_row]
            )
            right_remainder = -(
                characters[right_label, odd_column]
                + characters[right_label, even_column]
            )
            conditional = (
                left_at_row * right_at_column / n
                + left_remainder
                * right_remainder
                / (n * (n - 1))
            )
            by_descriptor = np.square(conditional / order).sum(axis=0)
            squared_sum += float(by_descriptor.sum())
            descriptor_counts[1] += len(grouped[1])
            maximum_per_descriptor[1] = max(
                maximum_per_descriptor[1],
                float(by_descriptor.max()),
            )

        if grouped[2]:
            data = np.asarray(grouped[2], dtype=np.int16)
            odd_row = data[:, 0][None, :]
            first_row = data[:, 1][None, :]
            first_label = data[:, 2][None, :]
            second_row = data[:, 3][None, :]
            second_label = data[:, 4][None, :]
            first_remainder = -(
                characters[first_label, odd_column]
                + characters[first_label, even_column]
            )
            second_remainder = -(
                characters[second_label, odd_column]
                + characters[second_label, even_column]
            )
            product_label = np.bitwise_xor(first_label, second_label)
            product_remainder = np.where(
                product_label == 0,
                order - 2,
                -(
                    characters[product_label, odd_column]
                    + characters[product_label, even_column]
                ),
            )
            right_pair_average = (
                first_remainder * second_remainder - product_remainder
            ) / ((n - 1) * (n - 2))
            left_remainder = -(
                characters[left_label, odd_row]
                + characters[left_label, first_row]
                + characters[left_label, second_row]
            )
            conditional = (
                characters[left_label, first_row]
                * characters[first_label, even_column]
                * second_remainder
                / (n * (n - 1))
                + characters[left_label, second_row]
                * characters[second_label, even_column]
                * first_remainder
                / (n * (n - 1))
                + left_remainder * right_pair_average / n
            )
            by_descriptor = np.square(conditional / order).sum(axis=0)
            squared_sum += float(by_descriptor.sum())
            descriptor_counts[2] += len(grouped[2])
            maximum_per_descriptor[2] = max(
                maximum_per_descriptor[2],
                float(by_descriptor.max()),
            )

    chained_coefficient = sqrt(
        squared_sum / (order**2 * (order - 1) ** 2)
    )
    return AdjacentJointSlice(
        order=order,
        fixed_singleton=fixed_singleton,
        fixed_triple=fixed_triple,
        cubic_extensions=len(cubic_supports),
        quintic_extensions=len(quintic_supports),
        descriptor_counts=tuple(descriptor_counts),
        maximum_per_descriptor=tuple(maximum_per_descriptor),
        squared_moment_sum=squared_sum,
        chained_coefficient=chained_coefficient,
    )


def joint_record_three_slice_sum(
    order: int,
    fixed_singleton: int,
    fixed_triple: tuple[int, int, int],
) -> AdjacentRecordThreeSlice:
    """Direct record-three M_35 slice, intended for q=4 and q=8."""

    if order > 8:
        raise ValueError("direct record-three slice is limited to q<=8")
    cubic_supports = supports_containing_with_records(
        order,
        (fixed_singleton,),
        3,
        row_record=1,
        column_record=3,
    )
    quintic_supports = supports_containing_with_records(
        order,
        fixed_triple,
        5,
        row_record=3,
    )
    squared_sum = 0.0
    for cubic in cubic_supports:
        squared_sum += sum(
            record_three_link_moment(order, cubic, quintic) ** 2
            for quintic in quintic_supports
        )
    chained_coefficient = sqrt(
        squared_sum / (order**2 * (order - 1) ** 2)
    )
    return AdjacentRecordThreeSlice(
        order=order,
        fixed_singleton=fixed_singleton,
        fixed_triple=fixed_triple,
        cubic_extensions=len(cubic_supports),
        quintic_extensions=len(quintic_supports),
        squared_moment_sum=squared_sum,
        chained_coefficient=chained_coefficient,
    )


def horizontal_record_three_slice_sum(
    order: int,
    *,
    chunk_size: int = 256,
) -> AdjacentRecordThreeSlice:
    """Vectorized record-three slice for a fixed horizontal triple."""

    if order > 16:
        raise ValueError("vectorized horizontal slice is limited to q<=16")
    fixed_singleton = 0
    fixed_triple = (0, 1, 2)
    cubic_supports = supports_containing_with_records(
        order,
        (fixed_singleton,),
        3,
        row_record=1,
        column_record=3,
    )
    quintic_supports = supports_containing_with_records(
        order,
        fixed_triple,
        5,
        row_record=3,
    )

    def odd_groups(
        support: tuple[int, ...],
        axis: int,
    ) -> tuple[tuple[int, ...], tuple[int, ...]]:
        groups: dict[int, list[int]] = {}
        for coordinate in support:
            labels = divmod(coordinate, order)
            groups.setdefault(labels[axis], []).append(labels[1 - axis])
        labels = tuple(
            sorted(
                label
                for label, values in groups.items()
                if len(values) % 2
            )
        )
        if len(labels) != 3:
            raise AssertionError(("record-three descriptor", support, axis))
        return labels, tuple(xor_values(groups[label]) for label in labels)

    left_groups = tuple(odd_groups(support, 1) for support in cubic_supports)
    right_groups = tuple(odd_groups(support, 0) for support in quintic_supports)
    left_labels = np.asarray([group[0] for group in left_groups], dtype=np.int16)
    left_xors = np.asarray([group[1] for group in left_groups], dtype=np.int16)
    right_labels = np.asarray([group[0] for group in right_groups], dtype=np.int16)
    right_xors = np.asarray([group[1] for group in right_groups], dtype=np.int16)
    characters = np.asarray(
        [
            [character(label, value) for value in range(order)]
            for label in range(order)
        ],
        dtype=np.int8,
    )
    squared_sum = 0.0
    denominator = order * (order - 1) * (order - 2)
    for start in range(0, len(quintic_supports), chunk_size):
        stop = min(start + chunk_size, len(quintic_supports))
        total = np.zeros((len(cubic_supports), stop - start), dtype=np.int16)
        for mapping in permutations(range(3)):
            phase = np.ones_like(total)
            for left_index, right_index in enumerate(mapping):
                phase *= characters[
                    left_xors[:, left_index, None],
                    right_labels[None, start:stop, right_index],
                ]
                phase *= characters[
                    left_labels[:, left_index, None],
                    right_xors[None, start:stop, right_index],
                ]
            total += phase
        squared_sum += float(np.square(total.astype(float)).sum())
    squared_sum /= denominator**2
    return AdjacentRecordThreeSlice(
        order=order,
        fixed_singleton=fixed_singleton,
        fixed_triple=fixed_triple,
        cubic_extensions=len(cubic_supports),
        quintic_extensions=len(quintic_supports),
        squared_moment_sum=squared_sum,
        chained_coefficient=sqrt(
            squared_sum / (order**2 * (order - 1) ** 2)
        ),
    )


def horizontal_adjacent_slice_certificate(
    order: int,
) -> HorizontalAdjacentSliceCertificate:
    """Closed horizontal-triple slices for both M_35 record sectors."""

    if order < 4 or order & (order - 1):
        raise ValueError(("power-of-two order at least four required", order))
    q = order
    no_even_group = comb(q - 3, 2) + q * (q - 1)
    one_even_group = (
        (q - 4) * q * (q - 1)
        + (q - 1) * comb(q, 2)
    )
    per_no_even = 3 / q**2
    per_one_even = 3 * (q**2 + 1) / (
        q**2 * (q - 1) ** 2
    )
    record_one = (
        no_even_group * per_no_even
        + one_even_group * per_one_even
    )
    # The record-three slice is D_1^star C_1:
    # [(q-1)(q-2)(3q-2)/2] * [3/((q-1)(q-2))].
    record_three = 3 * (3 * q - 2) / 2
    denominator = q**2 * (q - 1) ** 2
    return HorizontalAdjacentSliceCertificate(
        order=q,
        record_one_extensions=(no_even_group, one_even_group),
        record_one_squared_sum=float(record_one),
        record_three_squared_sum=float(record_three),
        record_one_coefficient=sqrt(record_one / denominator),
        record_three_coefficient=sqrt(record_three / denominator),
        combined_coefficient=sqrt(
            (record_one + record_three) / denominator
        ),
    )


def reduced_adjacent_orbit_witness(
    order: int,
    cubic_pair_difference: int,
    quintic_pair_difference: int,
    triple_shape: tuple[int, int, int],
    *,
    singleton: int = 0,
) -> ReducedAdjacentOrbitWitness:
    """Use M_11 orthogonality to diagonalize the first singleton link."""

    dimension = order * order
    cubic_pairs = pair_orbit(cubic_pair_difference, dimension)
    quintic_pairs = pair_orbit(quintic_pair_difference, dimension)
    triples = translated_sequence(triple_shape, dimension)
    row_count = len(cubic_pairs) * len(quintic_pairs)
    column_count = dimension * len(triples)
    row_indices: list[int] = []
    column_indices: list[int] = []
    values: list[float] = []

    for cubic_index, cubic_pair in enumerate(cubic_pairs):
        cubic_set = set(cubic_pair)
        for quintic_index, quintic_pair in enumerate(quintic_pairs):
            row = cubic_index * len(quintic_pairs) + quintic_index
            quintic_set = set(quintic_pair)
            for endpoint_singleton in range(dimension):
                if endpoint_singleton in cubic_set:
                    continue
                cubic_support = tuple(
                    sorted(cubic_pair + (endpoint_singleton,))
                )
                middle_value = record_one_link_moment(
                    order,
                    (singleton,),
                    cubic_support,
                )
                if middle_value == 0:
                    continue
                for triple_index, triple in enumerate(triples):
                    if quintic_set.intersection(triple):
                        continue
                    quintic_support = tuple(
                        sorted(quintic_pair + triple)
                    )
                    adjacent_value = record_one_link_moment(
                        order,
                        cubic_support,
                        quintic_support,
                    )
                    if adjacent_value == 0:
                        continue
                    row_indices.append(row)
                    column_indices.append(
                        endpoint_singleton * len(triples)
                        + triple_index
                    )
                    values.append(middle_value * adjacent_value)

    matrix = sparse.coo_matrix(
        (values, (row_indices, column_indices)),
        shape=(row_count, column_count),
    ).tocsr()
    gram = (matrix @ matrix.T).toarray()
    eigenvalues = np.linalg.eigvalsh((gram + gram.T) / 2)
    tolerance = 2e-11 * max(1.0, float(eigenvalues[-1]))
    if eigenvalues[0] < -tolerance:
        raise AssertionError(("non-PSD reduced row Gram", eigenvalues[0]))
    eigenvalues[eigenvalues < tolerance] = 0
    singular = np.sqrt(np.maximum(eigenvalues, 0))
    # M_11 is an N-by-N orthogonal matrix.  All singleton blocks are
    # translation-equivalent, so the full normalized coefficient is
    # N ||B||_* / (N^3/2) = 2 ||B||_*/N^2.
    return ReducedAdjacentOrbitWitness(
        order=order,
        cubic_pair_difference=cubic_pair_difference,
        quintic_pair_difference=quintic_pair_difference,
        triple_shape=triple_shape,
        block_rows=row_count,
        block_columns=column_count,
        nonzero=matrix.nnz,
        rank=int(np.count_nonzero(singular > 2e-10)),
        coefficient=float(2 * singular.sum() / dimension**2),
        normalized_operator=float(
            2 * singular[-1] / dimension**2
        ),
    )


def direct_adjacent_orbit_witness(
    order: int,
    cubic_pair_difference: int,
    quintic_pair_difference: int,
    triple_shape: tuple[int, int, int],
) -> AdjacentOrbitWitness:
    """Build and diagonalize one physical translation-orbit block."""

    moments = exact_link_moments(order)
    dimension = order * order
    cubic_pairs = pair_orbit(cubic_pair_difference, dimension)
    quintic_pairs = pair_orbit(quintic_pair_difference, dimension)
    triples = translated_sequence(triple_shape, dimension)
    three_index = {
        support: index
        for index, support in enumerate(moments.supports_three)
    }
    five_index = {
        support: index
        for index, support in enumerate(moments.supports_five)
    }

    row_count = dimension * len(cubic_pairs) * len(quintic_pairs)
    column_count = dimension * dimension * len(triples)
    row_indices: list[int] = []
    column_indices: list[int] = []
    values: list[float] = []

    for b in range(dimension):
        singleton_link = moments.moment_11[:, b]
        active_a = np.flatnonzero(singleton_link)
        for cubic_index, cubic_pair in enumerate(cubic_pairs):
            cubic_set = set(cubic_pair)
            for quintic_index, quintic_pair in enumerate(quintic_pairs):
                row = (
                    (b * len(cubic_pairs) + cubic_index)
                    * len(quintic_pairs)
                    + quintic_index
                )
                quintic_set = set(quintic_pair)
                for e in range(dimension):
                    if e in cubic_set:
                        continue
                    cubic_support = tuple(sorted(cubic_pair + (e,)))
                    cubic_support_index = three_index[cubic_support]
                    middle_value = moments.moment_13[
                        b,
                        cubic_support_index,
                    ]
                    if middle_value == 0:
                        continue
                    for triple_index, triple in enumerate(triples):
                        if quintic_set.intersection(triple):
                            continue
                        quintic_support = tuple(
                            sorted(quintic_pair + triple)
                        )
                        adjacent_value = moments.moment_35[
                            cubic_support_index,
                            five_index[quintic_support],
                        ]
                        if adjacent_value == 0:
                            continue
                        base_column = (
                            e * len(triples) + triple_index
                        )
                        for a in active_a:
                            column = (
                                a * dimension * len(triples)
                                + base_column
                            )
                            row_indices.append(row)
                            column_indices.append(column)
                            values.append(
                                singleton_link[a]
                                * middle_value
                                * adjacent_value
                            )

    matrix = sparse.coo_matrix(
        (values, (row_indices, column_indices)),
        shape=(row_count, column_count),
    ).tocsr()
    gram = (matrix @ matrix.T).toarray()
    eigenvalues = np.linalg.eigvalsh((gram + gram.T) / 2)
    tolerance = 2e-11 * max(1.0, float(eigenvalues[-1]))
    if eigenvalues[0] < -tolerance:
        raise AssertionError(("non-PSD row Gram", eigenvalues[0]))
    eigenvalues[eigenvalues < tolerance] = 0
    singular = np.sqrt(np.maximum(eigenvalues, 0))
    normalization = sqrt(row_count * column_count)
    return AdjacentOrbitWitness(
        order=order,
        cubic_pair_difference=cubic_pair_difference,
        quintic_pair_difference=quintic_pair_difference,
        triple_shape=triple_shape,
        rows=row_count,
        columns=column_count,
        nonzero=matrix.nnz,
        rank=int(np.count_nonzero(singular > 2e-10)),
        coefficient=float(singular.sum() / normalization),
        normalized_operator=float(singular[-1] / normalization),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        choices=("vertical", "l_shape"),
        default="vertical",
    )
    parser.add_argument(
        "--reduced-order",
        type=int,
        help="use the record-one formula and M_11 reduction",
    )
    parser.add_argument(
        "--joint-slice-order",
        type=int,
        help="evaluate the candidate record-aware singleton/triple slice",
    )
    parser.add_argument(
        "--record-three-slice-order",
        type=int,
        help="evaluate the direct record-three slice at q=4 or q=8",
    )
    parser.add_argument(
        "--horizontal-record-three-order",
        type=int,
        help="evaluate the vectorized horizontal record-three slice",
    )
    parser.add_argument(
        "--horizontal-slice-formula",
        type=int,
        help="evaluate the closed two-record horizontal slice",
    )
    arguments = parser.parse_args()
    if arguments.horizontal_slice_formula:
        result = horizontal_adjacent_slice_certificate(
            arguments.horizontal_slice_formula
        )
        print(
            "horizontal adjacent slice certificate: "
            f"q={result.order},"
            f"record_one_extensions={result.record_one_extensions},"
            f"record_one_squared_sum={result.record_one_squared_sum:.15g},"
            f"record_three_squared_sum="
            f"{result.record_three_squared_sum:.15g},"
            f"record_one_coefficient={result.record_one_coefficient:.15g},"
            f"record_three_coefficient="
            f"{result.record_three_coefficient:.15g},"
            f"combined_coefficient={result.combined_coefficient:.15g}"
        )
        return
    if arguments.horizontal_record_three_order:
        result = horizontal_record_three_slice_sum(
            arguments.horizontal_record_three_order
        )
        print(
            "horizontal adjacent record-three slice: "
            f"q={result.order},"
            f"cubic_extensions={result.cubic_extensions},"
            f"quintic_extensions={result.quintic_extensions},"
            f"squared_sum={result.squared_moment_sum:.15g},"
            f"chained_coefficient={result.chained_coefficient:.15g}"
        )
        return
    if arguments.record_three_slice_order:
        order = arguments.record_three_slice_order
        result = joint_record_three_slice_sum(
            order,
            0,
            (0, order, 2 * order),
        )
        print(
            "adjacent record-three M35 joint slice: "
            f"q={result.order},"
            f"singleton={result.fixed_singleton},"
            f"triple={result.fixed_triple},"
            f"cubic_extensions={result.cubic_extensions},"
            f"quintic_extensions={result.quintic_extensions},"
            f"squared_sum={result.squared_moment_sum:.15g},"
            f"chained_coefficient={result.chained_coefficient:.15g}"
        )
        return
    if arguments.joint_slice_order:
        order = arguments.joint_slice_order
        result = joint_record_one_slice_sum(
            order,
            0,
            (0, order, order + 1),
        )
        print(
            "adjacent M35 joint slice: "
            f"q={result.order},"
            f"singleton={result.fixed_singleton},"
            f"triple={result.fixed_triple},"
            f"cubic_extensions={result.cubic_extensions},"
            f"quintic_extensions={result.quintic_extensions},"
            f"descriptor_counts={result.descriptor_counts},"
            f"maximum_per_descriptor={result.maximum_per_descriptor},"
            f"squared_sum={result.squared_moment_sum:.15g},"
            f"chained_coefficient={result.chained_coefficient:.15g}"
        )
        return
    order = 4
    if arguments.case == "vertical":
        parameters = (order, order, (0, order, 2 * order))
    else:
        parameters = (order, order, (0, order, 1))
    if arguments.reduced_order:
        order = arguments.reduced_order
        if arguments.case == "vertical":
            parameters = (order, order, (0, order, 2 * order))
        else:
            parameters = (order, 1, (0, order, 1))
        reduced = reduced_adjacent_orbit_witness(order, *parameters)
        print(
            "reduced adjacent cubic-quintic orbit witness: "
            f"q={reduced.order},"
            f"x={reduced.cubic_pair_difference},"
            f"y={reduced.quintic_pair_difference},"
            f"triple={reduced.triple_shape},"
            f"block={reduced.block_rows}x{reduced.block_columns},"
            f"nonzero={reduced.nonzero},"
            f"rank={reduced.rank},"
            f"coefficient={reduced.coefficient:.15g},"
            f"normalized_operator={reduced.normalized_operator:.15g}"
        )
        return
    result = direct_adjacent_orbit_witness(order, *parameters)
    print(
        "adjacent cubic-quintic orbit witness: "
        f"q={result.order},"
        f"x={result.cubic_pair_difference},"
        f"y={result.quintic_pair_difference},"
        f"triple={result.triple_shape},"
        f"shape={result.rows}x{result.columns},"
        f"nonzero={result.nonzero},"
        f"rank={result.rank},"
        f"coefficient={result.coefficient:.15g},"
        f"normalized_operator={result.normalized_operator:.15g}"
    )


if __name__ == "__main__":
    main()
