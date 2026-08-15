#!/usr/bin/env python3
"""Small-q arbitrary-diagonal stress test for the alternating endpoint lift.

This is a computer-assisted q=2 result.  The exact orbit matrices and their
joint eigenspace multiplicities are checked over the integers.  A
supporting-hyperplane calculation for the jointly concave weighted trace
norm then bounds every pair of diagonal probability laws.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext

import numpy as np
import sympy as sp

from alternating_double_endpoint_spectrum import (
    pair_arrays,
    scaled_endpoint_weights,
)


@dataclass(frozen=True)
class WeightedQ2Certificate:
    orbit_mass: Decimal
    orbit_value: Decimal
    supporting_upper: float
    uniform_value: float
    high_columns: int
    low_columns: int


def sylvester(order: int) -> np.ndarray:
    result = np.array([[1.0]])
    while result.shape[0] < order:
        result = np.block([[result, result], [result, -result]])
    return result


def alternating_kernel_q2() -> tuple[np.ndarray, np.ndarray]:
    order = 2
    dimension = 4
    pair_left, pair_right = pair_arrays(dimension)
    pair_count = len(pair_left)
    hadamard = sylvester(dimension) / 2
    endpoint = np.empty((dimension, dimension, pair_count))
    for singleton in range(dimension):
        weights = scaled_endpoint_weights(
            order, singleton, pair_left, pair_right
        ).astype(float)
        endpoint[singleton] = (
            hadamard[singleton ^ pair_left ^ pair_right].T * weights
        )

    # Rows (i,b,d), columns (E,F,c).
    kernel = np.einsum(
        "ibe,bc,dcf->ibdecf",
        endpoint,
        hadamard,
        endpoint,
    ).transpose(0, 1, 2, 3, 5, 4).reshape(
        dimension**3, pair_count * pair_count * dimension
    )

    differences = pair_left ^ pair_right
    high_pair = np.empty((pair_count, pair_count), dtype=bool)
    for left_index, left_difference in enumerate(differences):
        for right_index, right_difference in enumerate(differences):
            high_pair[left_index, right_index] = (
                int(left_difference & right_difference).bit_count() % 2 == 0
            )
    high_columns = np.repeat(
        high_pair[:, :, None], dimension, axis=2
    ).reshape(-1)
    return kernel, high_columns


def check_exact_orbit_spectrum(
    kernel: np.ndarray, high_columns: np.ndarray
) -> None:
    """Verify the five joint high/low Gram eigenspaces exactly."""

    integer_kernel = np.rint(8 * kernel).astype(np.int64)
    if not np.array_equal(kernel, integer_kernel / 8):
        raise AssertionError("q=2 kernel is not an exact eighth-integer")

    high_gram = (
        integer_kernel[:, high_columns]
        @ integer_kernel[:, high_columns].T
    )
    low_gram = (
        integer_kernel[:, ~high_columns]
        @ integer_kernel[:, ~high_columns].T
    )
    if not np.array_equal(high_gram @ low_gram, low_gram @ high_gram):
        raise AssertionError("q=2 orbit Grams do not commute")

    identity = sp.eye(kernel.shape[0])
    exact_high = sp.Matrix(high_gram)
    exact_low = sp.Matrix(low_gram)
    expected_joint_spaces = (
        (0, 16, 24),
        (16, 0, 12),
        (16, 32, 12),
        (48, 0, 4),
        (16, 64, 12),
    )
    total = 0
    for high_eigenvalue, low_eigenvalue, multiplicity in expected_joint_spaces:
        stacked = sp.Matrix.vstack(
            exact_high - high_eigenvalue * identity,
            exact_low - low_eigenvalue * identity,
        )
        observed = kernel.shape[0] - stacked.rank()
        if observed != multiplicity:
            raise AssertionError(
                (
                    "joint orbit eigenspace",
                    high_eigenvalue,
                    low_eigenvalue,
                    observed,
                    multiplicity,
                )
            )
        total += observed
    if total != kernel.shape[0]:
        raise AssertionError(("joint orbit spectrum total", total))


def orbit_objective(mass: Decimal) -> Decimal:
    """Exact one-variable objective implied by the joint orbit spectrum."""

    one = Decimal(1)
    two = Decimal(2)
    with localcontext() as context:
        context.prec = 60
        root_three = Decimal(3).sqrt()
        return (
            Decimal(3)
            / (Decimal(8) * Decimal(6).sqrt())
            * (one - mass).sqrt()
            + (root_three + one) / Decimal(16) * mass.sqrt()
            + root_three / Decimal(16)
            + root_three / Decimal(16) * (two - mass).sqrt()
        )


def orbit_derivative(mass: Decimal) -> Decimal:
    one = Decimal(1)
    two = Decimal(2)
    with localcontext() as context:
        context.prec = 60
        root_three = Decimal(3).sqrt()
        return (
            -Decimal(3)
            / (
                Decimal(16)
                * Decimal(6).sqrt()
                * (one - mass).sqrt()
            )
            + (root_three + one)
            / (Decimal(32) * mass.sqrt())
            - root_three
            / (Decimal(32) * (two - mass).sqrt())
        )


def optimize_orbit_mass() -> tuple[Decimal, Decimal]:
    lower = Decimal("0.375")
    upper = Decimal("0.376")
    if orbit_derivative(lower) <= 0 or orbit_derivative(upper) >= 0:
        raise AssertionError("initial derivative bracket")
    for _ in range(180):
        middle = (lower + upper) / 2
        if orbit_derivative(middle) > 0:
            lower = middle
        else:
            upper = middle
    mass = (lower + upper) / 2
    return mass, orbit_objective(mass)


def supporting_hyperplane(
    kernel: np.ndarray,
    high_columns: np.ndarray,
    high_mass: float,
) -> tuple[float, float, float]:
    """Return value and a numerical global upper from joint concavity."""

    row_count, column_count = kernel.shape
    row_law = np.full(row_count, 1 / row_count)
    column_law = np.where(
        high_columns,
        high_mass / int(high_columns.sum()),
        (1 - high_mass) / int((~high_columns).sum()),
    )
    weighted = (
        np.sqrt(row_law)[:, None]
        * kernel
        * np.sqrt(column_law)[None, :]
    )
    left, singular_values, right = np.linalg.svd(
        weighted, full_matrices=False
    )
    polar = left @ right
    value = float(singular_values.sum())
    row_gradient = np.sum(polar * weighted, axis=1) / (2 * row_law)
    column_gradient = np.sum(polar * weighted, axis=0) / (
        2 * column_law
    )

    # For a jointly concave function, its tangent bounds every other pair
    # of simplex points.  Euler homogeneity cancels the base-point terms.
    upper = float(row_gradient.max() + column_gradient.max())
    gradient_spread = max(
        float(np.ptp(row_gradient)),
        float(np.ptp(column_gradient[high_columns])),
        float(np.ptp(column_gradient[~high_columns])),
        abs(
            float(column_gradient[high_columns].mean())
            - float(column_gradient[~high_columns].mean())
        ),
    )
    return value, upper, gradient_spread


def certificate() -> WeightedQ2Certificate:
    kernel, high_columns = alternating_kernel_q2()
    check_exact_orbit_spectrum(kernel, high_columns)
    mass, exact_value = optimize_orbit_mass()
    value, numerical_upper, gradient_spread = supporting_hyperplane(
        kernel, high_columns, float(mass)
    )
    if abs(Decimal(str(value)) - exact_value) > Decimal("3e-15"):
        raise AssertionError(("orbit objective", value, exact_value))
    if gradient_spread > 2e-15:
        raise AssertionError(("supporting gradient spread", gradient_spread))

    # Keep a conservative allowance far above observed double-precision
    # roundoff.  This is a certified numerical label, not interval arithmetic.
    supporting_upper = numerical_upper + 1e-12
    if supporting_upper >= 0.471845:
        raise AssertionError(("q=2 diagonal upper", supporting_upper))

    uniform_value = float(
        np.linalg.svd(kernel, compute_uv=False).sum()
        / np.sqrt(kernel.size)
    )
    return WeightedQ2Certificate(
        orbit_mass=mass,
        orbit_value=exact_value,
        supporting_upper=supporting_upper,
        uniform_value=uniform_value,
        high_columns=int(high_columns.sum()),
        low_columns=int((~high_columns).sum()),
    )


def main() -> None:
    result = certificate()
    print(
        "weighted alternating q=2 certificate passed: "
        f"high_columns={result.high_columns},"
        f"low_columns={result.low_columns},"
        f"high_orbit_mass={result.orbit_mass},"
        f"uniform={result.uniform_value:.15g},"
        f"optimum={result.orbit_value},"
        f"supporting_upper={result.supporting_upper:.15g}"
    )


if __name__ == "__main__":
    main()
