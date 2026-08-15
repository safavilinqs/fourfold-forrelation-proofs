#!/usr/bin/env python3
"""Regression for the transposed opposite-endpoint vertical witness."""

from __future__ import annotations

from itertools import combinations
from math import comb, sqrt
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from opposite_endpoint_orbit_scan import (
    cubic_weight,
    pair_orbit,
    quintic_weight,
    support_xor,
)
from opposite_endpoint_transposed_vertical_witness import (
    certified_dominant_class_witness,
    nystrom_vertical_witness,
    transposed_vertical_witness,
)


def character(left: int, right: int) -> int:
    return -1 if int(left & right).bit_count() % 2 else 1


def direct_q4() -> float:
    """Build one exact h-block before the symmetry quotient."""

    q = 4
    dimension = q * q
    pairs = tuple(
        pair
        for difference in range(q, dimension, q)
        for pair in pair_orbit(difference, dimension)
    )
    triples = tuple(
        tuple(row * q + column for row in rows)
        for column in range(q)
        for rows in combinations(range(q), 3)
    )
    block = np.zeros(
        (
            dimension * len(triples),
            len(pairs) * len(pairs),
        )
    )
    for endpoint_singleton in range(dimension):
        second_middle = endpoint_singleton
        for triple_index, triple in enumerate(triples):
            row = endpoint_singleton * len(triples) + triple_index
            for cubic_index, cubic_pair in enumerate(pairs):
                if endpoint_singleton in cubic_pair:
                    continue
                cubic = cubic_weight(
                    tuple(sorted(cubic_pair + (endpoint_singleton,))),
                    q,
                )
                for quintic_index, quintic_pair in enumerate(pairs):
                    if set(quintic_pair).intersection(triple):
                        continue
                    support = tuple(sorted(quintic_pair + triple))
                    quintic = quintic_weight(support, q)
                    phase = character(support_xor(support), second_middle)
                    column = cubic_index * len(pairs) + quintic_index
                    block[row, column] = (
                        cubic * quintic * phase / dimension
                    )
    block_nuclear = float(np.linalg.svd(block, compute_uv=False).sum())
    row_count = dimension**2 * len(triples)
    column_count = dimension * len(pairs) ** 2
    return dimension * block_nuclear / sqrt(row_count * column_count)


def main() -> None:
    q4 = transposed_vertical_witness(4)
    q8 = transposed_vertical_witness(8)
    q8_compressed = nystrom_vertical_witness(8)
    direct = direct_q4()
    if not np.isclose(q4.coefficient, direct, atol=3e-12):
        raise AssertionError(("direct q4 transposed witness", q4, direct))
    if not np.isclose(q4.coefficient, 0.0554262492783, atol=3e-12):
        raise AssertionError(q4)
    if not np.isclose(q8.coefficient, 0.0825349366921, atol=3e-12):
        raise AssertionError(q8)
    if q4.symmetry_block_size != 36 or q8.symmetry_block_size != 392:
        raise AssertionError(("symmetry block sizes", q4, q8))
    if not np.isclose(
        q8_compressed.coefficient_lower,
        0.0819670813997,
        atol=3e-12,
    ):
        raise AssertionError(("q8 compressed witness", q8_compressed))
    if q8_compressed.block_ranks != (24, 28, 24, 28, 24, 28, 24, 28):
        raise AssertionError(("q8 compressed ranks", q8_compressed))

    certified = certified_dominant_class_witness()
    if certified.retained_rank != 480 or certified.rational_trace_lower != 194:
        raise AssertionError(("q32 certificate dimensions", certified))
    if certified.contraction_row_upper >= 0.99981:
        raise AssertionError(("q32 contraction certificate", certified))
    if certified.spectral_error_upper >= 0.0015:
        raise AssertionError(("q32 spectral error", certified))
    if certified.computed_trace_lower <= 194.5:
        raise AssertionError(("q32 trace lower", certified))
    if not np.isclose(
        certified.coefficient_lower,
        0.0142810242047,
        atol=3e-13,
    ):
        raise AssertionError(("q32 coefficient lower", certified))
    print(
        "transposed opposite-endpoint witness passed: "
        f"q4={q4.coefficient:.12g},"
        f"q8={q8.coefficient:.12g},"
        f"q4_direct={direct:.12g},"
        f"q8_block_size={q8.symmetry_block_size},"
        f"q32_trace_lower={certified.rational_trace_lower},"
        f"q32_coefficient_lower={certified.coefficient_lower:.12g}"
    )


if __name__ == "__main__":
    main()
