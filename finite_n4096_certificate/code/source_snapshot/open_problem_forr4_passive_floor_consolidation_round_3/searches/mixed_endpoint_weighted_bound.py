#!/usr/bin/env python3
"""Arbitrary-diagonal bound for the worst mixed endpoint occurrence split.

The orientation has rows (i,b,F) and columns (E,c,d).  Translation twirling
reduces diagonal laws to distributions on the nonzero pair differences.
After endpoint Fourier transforms, every block is a diagonally weighted
Walsh matrix.  A three-type rank--Frobenius decomposition gives a closed
q-dependent upper bound without forming the huge physical matrix.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
from math import isqrt

import numpy as np

from alternating_double_endpoint_spectrum import (
    checked_fwht,
    orbit_correlation,
)


@dataclass(frozen=True)
class MixedEndpointBound:
    order: int
    amplitude_sums: tuple[Decimal, Decimal, Decimal]
    type_matrix: tuple[tuple[Decimal, ...], ...]
    row_sum_upper: Decimal
    spectral_upper: float
    old_fixed_split_bound: Decimal


@dataclass(frozen=True)
class RefinedMixedCertificate:
    vertical_vertical_upper: Fraction
    type_matrix_upper: tuple[tuple[Fraction, ...], ...]
    collatz_upper: Fraction
    advertised_upper: Fraction


def sylvester(order: int) -> np.ndarray:
    result = np.array([[1.0]])
    while result.shape[0] < order:
        result = np.block([[result, result], [result, -result]])
    return result


def mixed_fourier_value(
    order: int,
    row_difference_law: np.ndarray,
    column_difference_law: np.ndarray,
) -> float:
    """Exact reduced value for translation-invariant mixed-split laws."""

    dimension = order * order
    expected_shape = (dimension - 1,)
    row_law = np.asarray(row_difference_law, dtype=float)
    column_law = np.asarray(column_difference_law, dtype=float)
    if row_law.shape != expected_shape or column_law.shape != expected_shape:
        raise ValueError(("difference-law shape", expected_shape))
    if np.any(row_law < 0) or np.any(column_law < 0):
        raise ValueError("difference laws must be nonnegative")
    if not np.isclose(row_law.sum(), 1) or not np.isclose(
        column_law.sum(), 1
    ):
        raise ValueError("difference laws must have unit mass")

    endpoint_symbol = checked_fwht(
        orbit_correlation(order), axis=0
    )
    endpoint_fourier = np.sqrt(
        2 * endpoint_symbol[:, 1:].astype(float)
    ) / (order - 1)
    walsh = sylvester(dimension)[:, 1:]
    row_root = np.sqrt(row_law)
    column_root = np.sqrt(column_law)
    total = 0.0
    for row_frequency in range(dimension):
        row_weight = np.concatenate(
            ([0.0], row_root * endpoint_fourier[row_frequency])
        )
        for column_frequency in range(dimension):
            column_weight = (
                column_root * endpoint_fourier[column_frequency]
            )
            block = (
                row_weight[:, None]
                * walsh
                * column_weight[None, :]
            )
            total += float(
                np.linalg.svd(block, compute_uv=False).sum()
            )
    return total / dimension**3


def type_amplitude_sums(
    order: int,
) -> tuple[Decimal, Decimal, Decimal]:
    """Maximal summed endpoint Fourier RMS for V, H, and D pair types."""

    with localcontext() as context:
        context.prec = 70
        q = Decimal(order)
        n = q - 1
        endpoint_large = (q * q - 2 * q + 2) / n
        half_space = (2 * (q - 2) / n).sqrt()
        diagonal_half_space = (
            (2 * (q * q - 2 * q + 2)).sqrt() / (n * n)
        )
        vertical = (
            2
            + n * endpoint_large
            + n * half_space
            + n * n * half_space
        )
        horizontal = (
            2
            + n * half_space
            + n * (2 / n)
            + n * n * (half_space / n)
        )
        diagonal = (
            2
            + n * half_space
            + n * (half_space / n)
            + n * n * diagonal_half_space
        )
        return vertical, horizontal, diagonal


def sqrt_fraction_upper(
    value: Fraction, decimal_places: int = 60
) -> Fraction:
    """Exact rational upper approximation to a nonnegative square root."""

    if value < 0:
        raise ValueError(value)
    scale = 10**decimal_places
    target = value.numerator * scale * scale
    root = isqrt(target // value.denominator)
    if root * root * value.denominator < target:
        root += 1
    result = Fraction(root, scale)
    if result * result < value:
        raise AssertionError(("square-root upper", value, result))
    return result


def vertical_vertical_coefficient(order: int) -> Decimal:
    """Exact-radical value for uniform laws on the two vertical types."""

    if order < 4 or order & (order - 1):
        raise ValueError(("order must be a power of two at least four", order))
    with localcontext() as context:
        context.prec = 70
        q = Decimal(order)
        n = q - 1
        half = q / 2
        hyperplane = half - 1
        constant_sum = q * q - 2 * q + 4
        root_q = q.sqrt()
        full_nuclear = 1 + (q - 2) * root_q
        constant_hyperplane = (
            (half - 2) * root_q + (half + 1).sqrt()
        )
        dot_zero = (
            (q / 4 - 2) * root_q + (2 * q + 1).sqrt()
        )
        dot_one = 1 + (half - 2) * half.sqrt()
        numerator = (
            constant_sum * constant_sum / n * full_nuclear
            + 4 * q * constant_sum * constant_hyperplane
            + 4
            * q
            * q
            * (hyperplane * dot_zero + half * dot_one)
        )
        return numerator / q**6


def vertical_vertical_upper_q32() -> Fraction:
    """Exact-rational upper for the vertical/vertical type coefficient."""

    q = 32
    n = q - 1
    half = q // 2
    hyperplane = half - 1
    constant_sum = q * q - 2 * q + 4
    root_q = sqrt_fraction_upper(Fraction(q))
    full_nuclear = 1 + (q - 2) * root_q
    constant_hyperplane = (
        (half - 2) * root_q
        + sqrt_fraction_upper(Fraction(half + 1))
    )
    dot_zero = (
        (q // 4 - 2) * root_q
        + sqrt_fraction_upper(Fraction(2 * q + 1))
    )
    dot_one = (
        1
        + (half - 2) * sqrt_fraction_upper(Fraction(half))
    )
    numerator = (
        Fraction(constant_sum * constant_sum, n) * full_nuclear
        + 4 * q * constant_sum * constant_hyperplane
        + 4
        * q
        * q
        * (hyperplane * dot_zero + half * dot_one)
    )
    return numerator / q**6


def refined_q32_certificate() -> RefinedMixedCertificate:
    """Refine the three-type bound using the exact V/V contraction.

    Every radical is replaced by an exact rational upper.  Collatz--Wielandt
    with an explicit positive rational vector then certifies the Perron
    root of the entrywise upper matrix.
    """

    q = 32
    n = q - 1
    dimension = q * q
    endpoint_large = Fraction(q * q - 2 * q + 2, n)
    half_space = sqrt_fraction_upper(Fraction(2 * (q - 2), n))
    diagonal_half_space = (
        sqrt_fraction_upper(Fraction(2 * (q * q - 2 * q + 2)))
        / (n * n)
    )
    amplitudes = (
        2 + n * endpoint_large + n * half_space + n * n * half_space,
        4 + 2 * n * half_space,
        2 + n * half_space + half_space + n * n * diagonal_half_space,
    )
    counts = (n, n, n * n)
    matrix = [
        [
            sqrt_fraction_upper(Fraction(min(counts[row], counts[column])))
            * amplitudes[row]
            * amplitudes[column]
            / dimension**3
            for column in range(3)
        ]
        for row in range(3)
    ]
    vertical_vertical = vertical_vertical_upper_q32()
    matrix[0][0] = vertical_vertical

    vector = (
        Fraction(1),
        Fraction(54155, 10**6),
        Fraction(54752, 10**6),
    )
    ratios = tuple(
        sum(matrix[row][column] * vector[column] for column in range(3))
        / vector[row]
        for row in range(3)
    )
    collatz_upper = max(ratios)
    advertised_upper = Fraction(20343, 10**6)
    if collatz_upper >= advertised_upper:
        raise AssertionError(
            ("refined mixed Collatz certificate", collatz_upper)
        )
    return RefinedMixedCertificate(
        vertical_vertical_upper=vertical_vertical,
        type_matrix_upper=tuple(tuple(row) for row in matrix),
        collatz_upper=collatz_upper,
        advertised_upper=advertised_upper,
    )


def bound(order: int) -> MixedEndpointBound:
    with localcontext() as context:
        context.prec = 70
        q = Decimal(order)
        dimension = q * q
        counts = (q - 1, q - 1, (q - 1) ** 2)
        amplitudes = type_amplitude_sums(order)
        matrix = tuple(
            tuple(
                min(counts[row], counts[column]).sqrt()
                * amplitudes[row]
                * amplitudes[column]
                / dimension**3
                for column in range(3)
            )
            for row in range(3)
        )
        row_sums = tuple(sum(row) for row in matrix)
        row_sum_upper = max(row_sums)
        energy_one = (q * q + 2) / (2 * q * q)
        energy_two = (q * q - 2 * q + 2) / (
            q * q * (q - 1)
        )
        old_bound = (energy_one * energy_two).sqrt()

    float_matrix = np.array(
        [[float(entry) for entry in row] for row in matrix]
    )
    spectral_upper = float(np.linalg.norm(float_matrix, 2) + 1e-15)
    if spectral_upper > float(row_sum_upper) + 2e-15:
        raise AssertionError(("matrix norm versus row sum", order))
    return MixedEndpointBound(
        order=order,
        amplitude_sums=amplitudes,
        type_matrix=matrix,
        row_sum_upper=row_sum_upper,
        spectral_upper=spectral_upper,
        old_fixed_split_bound=old_bound,
    )


def improved_deterministic_ledger() -> tuple[float, float, float]:
    """Safe local-triangle ledger for occupation (2,1,1,2) at q=32.

    The same-middle masks now use their exact arbitrary-diagonal
    coefficients.  The alternating masks retain this file's analytic
    mixed-orientation upper and the companion same-orientation certificate.
    """

    from same_middle_weighted_bound import deterministic_ledger

    refined = refined_q32_certificate().advertised_upper
    raw, attenuated, remaining_margin = deterministic_ledger(
        alternating_equal=Decimal("0.010905"),
        alternating_mixed=(
            Decimal(refined.numerator) / Decimal(refined.denominator)
        ),
    )
    return float(raw), float(attenuated), float(remaining_margin)


def main() -> None:
    for order in (4, 8, 16, 32):
        result = bound(order)
        ratio = result.row_sum_upper / result.old_fixed_split_bound
        print(
            f"q={order},N={order * order},"
            f"row_sum_upper={result.row_sum_upper},"
            f"spectral_upper={result.spectral_upper:.15g},"
            f"old_bound={result.old_fixed_split_bound},"
            f"improvement_ratio={ratio}"
        )
    refined = refined_q32_certificate()
    print(
        "q=32 refined mixed certificate: "
        f"vertical_vertical_upper={float(refined.vertical_vertical_upper):.15g},"
        f"collatz_upper={float(refined.collatz_upper):.15g},"
        f"advertised_upper={float(refined.advertised_upper):.15g}"
    )
    raw, attenuated, margin = improved_deterministic_ledger()
    print(
        "same-middle-refined fixed-split triangle ledger: "
        f"raw={raw:.15g},attenuated={attenuated:.15g},"
        f"available_margin={margin:.15g},"
        f"slack={margin-attenuated:.15g}"
    )


if __name__ == "__main__":
    main()
