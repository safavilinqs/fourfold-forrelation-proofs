#!/usr/bin/env python3
"""Regression for the exact two-complete-flag collective obstruction."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from active_two_flag_collective_obstruction import (  # noqa: E402
    EXPECTED_SQUARED_SPECTRUM,
    chain_numerator,
    endpoint_factor_data,
    fixed_unit_factor,
    sylvester_sign,
    two_flag_collective_audit,
)


def folded_states(blocks: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    dimension = blocks.shape[1]
    hadamard = sylvester_sign(dimension) / np.sqrt(dimension)
    uniform = np.ones(dimension) / np.sqrt(dimension)
    left = blocks[1] * (hadamard @ (blocks[0] * uniform))
    right = hadamard @ (
        blocks[2] * (hadamard @ (blocks[3] * uniform))
    )
    return left, right


def tensor_product_checks() -> None:
    """Protect endpoint multiplicativity and the isometric N=1024 lift."""

    data = endpoint_factor_data()
    hadamard_16 = sylvester_sign(16)
    samples = (
        (1, 1, 0, 1),
        (1, -1, 7, 11),
        (-1, 1, 13, 17),
        (-1, -1, 29, 31),
    )
    for first_sign, second_sign, first_index, second_index in samples:
        first = data[first_sign][0][first_index]
        second = data[second_sign][0][second_index]
        tensor_blocks = np.asarray(
            [
                np.kron(first[block], second[block])
                for block in range(4)
            ]
        )
        numerator = chain_numerator(tensor_blocks, hadamard_16)
        expected_value = first_sign * second_sign / 4
        actual_value = numerator / (16**2.5)
        if not np.isclose(actual_value, expected_value, atol=2e-14):
            raise AssertionError(
                ("tensor endpoint multiplicativity", actual_value)
            )

        left_first, right_first = folded_states(first)
        left_second, right_second = folded_states(second)
        left_tensor, right_tensor = folded_states(tensor_blocks)
        if not np.allclose(
            left_tensor,
            np.kron(left_first, left_second),
            atol=2e-14,
        ):
            raise AssertionError("left folded state is not multiplicative")
        if not np.allclose(
            right_tensor,
            np.kron(right_first, right_second),
            atol=2e-14,
        ):
            raise AssertionError("right folded state is not multiplicative")

    unit = fixed_unit_factor()
    left_unit, right_unit = folded_states(unit)
    if not np.isclose(left_unit @ right_unit, 1, atol=2e-14):
        raise AssertionError("fixed lift factor does not have unit overlap")
    if not np.allclose(left_unit, right_unit, atol=2e-14):
        raise AssertionError("unit-overlap lift is not a common tensor factor")


def main() -> None:
    result = two_flag_collective_audit()
    if result.squared_singular_spectrum != EXPECTED_SQUARED_SPECTRUM:
        raise AssertionError(("squared singular spectrum", result))
    expected_distance = (
        sp.Rational(7, 88)
        + 15 * sp.sqrt(5) / 242
        + 9 * sp.sqrt(41) / 968
    )
    if sp.simplify(result.trace_distance_exact - expected_distance) != 0:
        raise AssertionError(("exact trace distance", result))
    if not result.trace_distance < 1 / 3:
        raise AssertionError(("two flags unexpectedly distinguish", result))
    if not result.helstrom_error > 1 / 3:
        raise AssertionError(("two-flag error threshold", result))
    if result.nonzero_rank != 100:
        raise AssertionError(("two-flag rank", result))
    if (result.gram_components, result.gram_component_size) != (16, 32):
        raise AssertionError(("Gram decomposition", result))
    tensor_product_checks()
    print(
        "active two-flag collective obstruction passed: "
        f"N0={result.tensor_endpoint_dimension},"
        f"N={result.target_dimension},"
        f"factor_inputs={result.endpoint_inputs_per_sign},"
        f"rank={result.nonzero_rank},"
        f"trace_distance={result.trace_distance:.12g},"
        f"helstrom_error={result.helstrom_error:.12g},"
        f"error_margin={result.error_margin_over_one_third:.12g}"
    )


if __name__ == "__main__":
    main()
