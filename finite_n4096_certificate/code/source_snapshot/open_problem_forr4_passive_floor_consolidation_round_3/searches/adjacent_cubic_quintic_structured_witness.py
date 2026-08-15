#!/usr/bin/env python3
"""Structured adjacent mixed-orbit witness through moderate q.

The row law is uniform on

    cubic selected-pair XOR x = (0,c), c != 0,
    quintic selected-pair XOR y = (r,0), r != 0,

and the column law is the translation orbit of one vertical triple.  This is
the largest simple subfamily visible in the exact q=4 mixed-law optimizer.
Because the column shape is fixed, every frequency block has one column and
its nuclear norm is its Euclidean norm.  The script accumulates these norms
without storing the full mixed-orbit tensor.

The computation uses exact signed-permutation link formulas followed by fast
Walsh transforms.  Values beyond q=4 are numerical physical lower witnesses,
not arbitrary-law upper bounds or scaling theorems.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np

from adjacent_cubic_quintic_mixed_orbit_q4 import combined_link_moment
from adjacent_cubic_quintic_orbit_witness import (
    parity_record_size,
    record_one_link_moment,
    xor_values,
)


@dataclass(frozen=True)
class StructuredAdjacentWitness:
    order: int
    row_types: int
    triple: tuple[int, int, int]
    coefficient: float
    record_one_coefficient: float
    record_three_coefficient: float
    interference_gain: float


def fwht_axis(values: np.ndarray, axis: int) -> np.ndarray:
    """Unnormalized Walsh transform along one power-of-two axis."""

    result = np.asarray(values, dtype=float).copy()
    result = np.swapaxes(result, axis, -1)
    length = result.shape[-1]
    if length == 0 or length & (length - 1):
        raise ValueError(("Walsh axis length", length))
    width = 1
    while width < length:
        shaped = result.reshape(*result.shape[:-1], -1, 2, width)
        left = shaped[..., 0, :].copy()
        right = shaped[..., 1, :].copy()
        shaped[..., 0, :] = left + right
        shaped[..., 1, :] = left - right
        width *= 2
    return np.swapaxes(result, axis, -1)


def twisted_walsh(kernel: np.ndarray) -> np.ndarray:
    dimension = kernel.shape[0]
    characters = np.empty_like(kernel)
    for left in range(dimension):
        for right in range(dimension):
            characters[left, right] = (
                -1 if int(left & right).bit_count() % 2 else 1
            )
    return fwht_axis(fwht_axis(kernel * characters, 0), 1)


def vertical_triple(order: int) -> tuple[int, int, int]:
    if order < 4:
        raise ValueError(order)
    return (0, order, 2 * order)


def structured_adjacent_witness(order: int) -> StructuredAdjacentWitness:
    if order < 4 or order & (order - 1):
        raise ValueError(("power-of-two order at least four required", order))
    dimension = order**2
    triple = vertical_triple(order)
    triple_xor = xor_values(list(triple))
    horizontal = tuple(range(1, order))
    vertical = tuple(row * order for row in range(1, order))
    squared_total = np.zeros((dimension, dimension), dtype=float)
    squared_record_one = np.zeros_like(squared_total)
    squared_record_three = np.zeros_like(squared_total)

    for x in horizontal:
        middle = np.zeros(dimension, dtype=float)
        cubics: list[tuple[int, int, int] | None] = []
        cubic_records = np.zeros(dimension, dtype=np.int8)
        for s in range(dimension):
            if s in (0, x):
                cubics.append(None)
                continue
            cubic = tuple(sorted((0, x, s)))
            cubics.append(cubic)
            middle[s] = record_one_link_moment(order, (0,), cubic)
            cubic_records[s] = parity_record_size(order, cubic, axis=1)
        for y in vertical:
            record_one = np.zeros((dimension, dimension), dtype=float)
            record_three = np.zeros_like(record_one)
            for s, cubic in enumerate(cubics):
                if cubic is None or middle[s] == 0:
                    continue
                for t in range(dimension):
                    shifted = tuple(value ^ t for value in triple)
                    if 0 in shifted or y in shifted:
                        continue
                    quintic = tuple(sorted((0, y) + shifted))
                    adjacent = combined_link_moment(
                        order,
                        cubic,
                        quintic,
                    )
                    value = middle[s] * adjacent
                    if cubic_records[s] == 1:
                        record_one[s, t] = value
                    elif cubic_records[s] == 3:
                        record_three[s, t] = value
            spectrum_one = twisted_walsh(record_one)
            spectrum_three = twisted_walsh(record_three)
            spectrum = spectrum_one + spectrum_three
            shifted_indices = np.bitwise_xor(
                np.arange(dimension),
                triple_xor,
            )
            squared_total += np.square(spectrum[shifted_indices])
            squared_record_one += np.square(spectrum_one[shifted_indices])
            squared_record_three += np.square(spectrum_three[shifted_indices])

    scale = 1 / (dimension**2 * (order - 1))
    coefficient = scale * float(np.sqrt(squared_total).sum())
    coefficient_one = scale * float(np.sqrt(squared_record_one).sum())
    coefficient_three = scale * float(np.sqrt(squared_record_three).sum())
    separate = coefficient_one + coefficient_three
    return StructuredAdjacentWitness(
        order=order,
        row_types=(order - 1) ** 2,
        triple=triple,
        coefficient=coefficient,
        record_one_coefficient=coefficient_one,
        record_three_coefficient=coefficient_three,
        interference_gain=(coefficient / separate if separate else 0.0),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orders", nargs="*", type=int, default=(4, 8))
    arguments = parser.parse_args()
    for order in arguments.orders:
        result = structured_adjacent_witness(order)
        print(
            "structured adjacent mixed witness: "
            f"q={result.order},N={result.order**2},"
            f"row_types={result.row_types},"
            f"triple={result.triple},"
            f"record_one={result.record_one_coefficient:.15g},"
            f"record_three={result.record_three_coefficient:.15g},"
            f"combined={result.coefficient:.15g},"
            f"combined_over_separate={result.interference_gain:.12g}"
        )


if __name__ == "__main__":
    main()
