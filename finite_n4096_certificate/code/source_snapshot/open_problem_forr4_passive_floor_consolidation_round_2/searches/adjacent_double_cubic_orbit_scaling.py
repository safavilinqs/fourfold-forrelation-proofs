#!/usr/bin/env python3
"""Exact q=4,8 scaling of the worst adjacent double-cubic slice."""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from math import comb, sqrt

import numpy as np
from scipy import sparse


def endpoint_supports(order: int):
    supports = []
    type_two = []
    alpha = []
    delta_column = []
    for odd_column in range(order):
        for rows in combinations(range(order), 3):
            supports.append(tuple(sorted(row * order + odd_column for row in rows)))
            type_two.append(False)
            alpha.append(0)
            delta_column.append(0)
    for odd_column in range(order):
        for even_column in range(order):
            if even_column == odd_column:
                continue
            for odd_row in range(order):
                for even_rows in combinations(range(order), 2):
                    supports.append(
                        tuple(
                            sorted(
                                (odd_row * order + odd_column,)
                                + tuple(
                                    row * order + even_column
                                    for row in even_rows
                                )
                            )
                        )
                    )
                    type_two.append(True)
                    alpha.append(even_rows[0] ^ even_rows[1])
                    delta_column.append(odd_column ^ even_column)
    if len(set(supports)) != len(supports):
        raise AssertionError("duplicate endpoint supports")
    return (
        supports,
        np.array(type_two, dtype=bool),
        np.array(alpha, dtype=np.int16),
        np.array(delta_column, dtype=np.int16),
    )


def l_supports(order: int):
    supports = []
    delta_row = []
    beta = []
    for odd_row in range(order):
        for even_row in range(order):
            if even_row == odd_row:
                continue
            for odd_column in range(order):
                for even_column in range(order):
                    if even_column == odd_column:
                        continue
                    supports.append(
                        tuple(
                            sorted(
                                (
                                    odd_row * order + even_column,
                                    even_row * order + odd_column,
                                    even_row * order + even_column,
                                )
                            )
                        )
                    )
                    delta_row.append(odd_row ^ even_row)
                    beta.append(odd_column ^ even_column)
    if len(set(supports)) != len(supports):
        raise AssertionError("duplicate L supports")
    return (
        supports,
        np.array(delta_row, dtype=np.int16),
        np.array(beta, dtype=np.int16),
    )


def part_lists(supports, size: int):
    result = defaultdict(list)
    for index, support in enumerate(supports):
        for part in combinations(support, size):
            result[part].append(index)
    keys = sorted(result)
    return keys, [np.array(result[key], dtype=np.int32) for key in keys]


def sparse_part_incidence(lists, support_count: int):
    row = []
    column = []
    for row_index, values in enumerate(lists):
        row.extend([row_index] * len(values))
        column.extend(values.tolist())
    data = np.ones(len(row), dtype=np.int16)
    return sparse.csr_matrix(
        (data, (row, column)),
        shape=(len(lists), support_count),
    )


def exact_worst(order: int):
    endpoint, type_two, alpha, delta_column = endpoint_supports(order)
    middle, delta_row, beta = l_supports(order)
    orthogonal = np.empty((order, order), dtype=bool)
    for left in range(order):
        for right in range(order):
            orthogonal[left, right] = (left & right).bit_count() % 2 == 0
    exceptional = (
        type_two[:, None]
        & orthogonal[alpha[:, None], delta_row[None, :]]
        & orthogonal[beta[None, :], delta_column[:, None]]
    )

    baseline = 1 / (order * order * (order - 1) ** 2)
    enhanced = (order + 2) ** 2 / (
        order * order * (order - 1) ** 2 * (order - 2) ** 2
    )
    increment = enhanced - baseline

    endpoint_parts = {}
    middle_parts = {}
    middle_incidence = {}
    middle_sizes = {}
    for size in range(4):
        endpoint_parts[size] = part_lists(endpoint, size)
        middle_parts[size] = part_lists(middle, size)
        middle_incidence[size] = sparse_part_incidence(
            middle_parts[size][1], len(middle)
        )
        middle_sizes[size] = np.array(
            [len(values) for values in middle_parts[size][1]]
        )

    slice_maximum = {}
    slice_key = {}
    for endpoint_size in range(4):
        endpoint_keys, endpoint_lists = endpoint_parts[endpoint_size]
        for middle_size in range(4):
            middle_keys, _ = middle_parts[middle_size]
            maximum = 0.0
            maximum_key = None
            for key, endpoint_indices in zip(endpoint_keys, endpoint_lists):
                exceptional_by_middle = exceptional[endpoint_indices].sum(axis=0)
                exceptional_counts = (
                    middle_incidence[middle_size] @ exceptional_by_middle
                )
                values = (
                    baseline
                    * len(endpoint_indices)
                    * middle_sizes[middle_size]
                    + increment * exceptional_counts
                )
                location = int(np.argmax(values))
                if values[location] > maximum:
                    maximum = float(values[location])
                    maximum_key = (key, middle_keys[location])
            slice_maximum[(endpoint_size, middle_size)] = maximum
            slice_key[(endpoint_size, middle_size)] = maximum_key

    coefficients = {}
    worst = (0.0, None)
    for endpoint_size in range(4):
        for middle_size in range(4):
            row_energy = (
                slice_maximum[(endpoint_size, middle_size)]
                / (order - 1) ** 2
            )
            column_energy = (
                slice_maximum[(3 - endpoint_size, 3 - middle_size)]
                / (order * order * (order - 1) ** 2)
            )
            coefficient = min(np.sqrt(row_energy), np.sqrt(column_energy))
            coefficients[(endpoint_size, middle_size)] = (
                coefficient,
                row_energy,
                column_energy,
            )
            if coefficient > worst[0]:
                worst = (coefficient, (endpoint_size, middle_size))

    maximum_f = slice_maximum[(2, 2)]
    maximum_g = slice_maximum[(1, 1)]
    row_energy = coefficients[(2, 2)][1]
    column_energy = coefficients[(2, 2)][2]
    coefficient = coefficients[(2, 2)][0]
    return {
        "order": order,
        "endpoint_supports": len(endpoint),
        "middle_supports": len(middle),
        "F22": maximum_f,
        "G11": maximum_g,
        "row_energy": row_energy,
        "column_energy": column_energy,
        "coefficient": coefficient,
        "F_key": slice_key[(2, 2)],
        "G_key": slice_key[(1, 1)],
        "coefficients": coefficients,
        "worst": worst,
    }


def exact_two_two_formulas(order: int):
    """Closed counts for the two slice maxima used by the 2|1 split."""
    f22 = (
        2 * (order * order - 2) / (order * order * (order - 1))
        + 8 / ((order - 1) * (order - 2))
    )
    g11 = 9 + (
        3
        * (order - 1)
        * (order - 2 + 3 * order * (order - 1))
        / (2 * order * order)
    )
    return f22, g11


def uniform_incidence_bounds(order: int):
    """All-q fixed-split bounds using the largest squared link entry."""
    endpoint_degrees = [
        order * comb(order, 3)
        + order * order * (order - 1) * comb(order, 2),
        comb(order - 1, 2)
        + (order - 1) * comb(order, 2)
        + order * (order - 1) ** 2,
        order * order - 2,
        1,
    ]
    middle_degrees = [
        order * order * (order - 1) ** 2,
        3 * (order - 1) ** 2,
        2 * (order - 1),
        1,
    ]
    maximum_squared_entry = (order + 2) ** 2 / (
        order * order * (order - 1) ** 2 * (order - 2) ** 2
    )
    coefficients = {}
    worst = (0.0, None)
    for endpoint_size in range(4):
        for middle_size in range(4):
            row = sqrt(
                maximum_squared_entry
                * endpoint_degrees[endpoint_size]
                * middle_degrees[middle_size]
            ) / (order - 1)
            column = sqrt(
                maximum_squared_entry
                * endpoint_degrees[3 - endpoint_size]
                * middle_degrees[3 - middle_size]
            ) / (order * (order - 1))
            coefficient = min(row, column)
            coefficients[(endpoint_size, middle_size)] = coefficient
            if coefficient > worst[0]:
                worst = (coefficient, (endpoint_size, middle_size))
    return coefficients, worst


def main() -> None:
    rows = [exact_worst(order) for order in (4, 8)]
    expected_q4 = np.sqrt(35 / 256)
    if not np.isclose(rows[0]["coefficient"], expected_q4, atol=2e-12):
        raise AssertionError(("q=4 adjacent slice", rows[0], expected_q4))
    for row in rows:
        expected_f, expected_g = exact_two_two_formulas(row["order"])
        if not np.isclose(row["F22"], expected_f, atol=2e-12):
            raise AssertionError(("F22 formula", row, expected_f))
        if not np.isclose(row["G11"], expected_g, atol=2e-12):
            raise AssertionError(("G11 formula", row, expected_g))
        if row["worst"][1] != (2, 2):
            raise AssertionError(("finite split maximum", row))
    q32_coefficients, q32_worst = uniform_incidence_bounds(32)
    if q32_worst[1] != (2, 2):
        raise AssertionError(("q=32 uniform split maximum", q32_worst))
    print("adjacent double-cubic orbit scaling:")
    for row in rows:
        print(
            f"q={row['order']},endpoint_supports={row['endpoint_supports']},"
            f"middle_supports={row['middle_supports']},F22={row['F22']:.12g},"
            f"G11={row['G11']:.12g},row_energy={row['row_energy']:.12g},"
            f"column_energy={row['column_energy']:.12g},"
            f"coefficient={row['coefficient']:.12g},"
            f"F_key={row['F_key']},G_key={row['G_key']},"
            f"worst={row['worst']}"
        )
        for placement, values in sorted(row["coefficients"].items()):
            print(
                f"  split={placement},coefficient={values[0]:.12g},"
                f"row={values[1]:.12g},column={values[2]:.12g}"
            )
    f32, g32 = exact_two_two_formulas(32)
    exact_two_two = min(sqrt(f32) / 31, sqrt(g32) / (32 * 31))
    print(
        "q=32 uniform incidence bounds: "
        f"worst={q32_worst[0]:.12g},placement={q32_worst[1]},"
        f"exact_2_2={exact_two_two:.12g},"
        f"attenuated_worst={q32_worst[0] * (5 / 6) ** 8:.12g}"
    )
    for placement, coefficient in sorted(q32_coefficients.items()):
        print(f"  split={placement},bound={coefficient:.12g}")


if __name__ == "__main__":
    main()
