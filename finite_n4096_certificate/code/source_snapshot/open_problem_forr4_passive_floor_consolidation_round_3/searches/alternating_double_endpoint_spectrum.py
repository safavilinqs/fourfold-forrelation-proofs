#!/usr/bin/env python3
"""Exact Walsh certificate for the alternating double-endpoint row Gram.

The cubic endpoint coefficient has the form

    A(i; {u,v}, b) = w_i({u,v}) H_N(i xor u xor v, b),

where (q-1) w is integer-valued.  After a Walsh transform in the first
middle label, every row-Gram block is a two-dimensional XOR convolution.
This script constructs its complete spectrum using integer orbit counts and
integer fast Walsh transforms.  It never forms the N^3 by N^3 row Gram.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext

import numpy as np


INT64_MAX = np.iinfo(np.int64).max


@dataclass(frozen=True)
class SpectrumCertificate:
    order: int
    dimension: int
    denominator: int
    numerators: np.ndarray
    multiplicities: np.ndarray
    coefficient: Decimal


def checked_fwht(values: np.ndarray, axis: int) -> np.ndarray:
    """Unnormalized Walsh transform, with an int64 overflow guard."""

    transformed = np.array(values, dtype=np.int64, copy=True)
    transformed = np.swapaxes(transformed, axis, -1)
    size = transformed.shape[-1]
    if size < 1 or size & (size - 1):
        raise ValueError(("Walsh axis must have power-of-two length", size))

    width = 1
    while width < size:
        maximum = int(np.max(np.abs(transformed)))
        if maximum > INT64_MAX // 2:
            raise OverflowError(("int64 Walsh stage", width, maximum))
        blocks = transformed.reshape(
            *transformed.shape[:-1], -1, 2 * width
        )
        left = blocks[..., :width].copy()
        right = blocks[..., width:].copy()
        blocks[..., :width] = left + right
        blocks[..., width:] = left - right
        width *= 2
    return np.swapaxes(transformed, axis, -1)


def pair_arrays(dimension: int) -> tuple[np.ndarray, np.ndarray]:
    left, right = np.triu_indices(dimension, k=1)
    return left.astype(np.int32), right.astype(np.int32)


def scaled_endpoint_weights(
    order: int,
    singleton: int,
    pair_left: np.ndarray,
    pair_right: np.ndarray,
) -> np.ndarray:
    """Return (q-1) w_i(E), exactly in {q-1,-1,0}.

    The value is q-1 when all three endpoint coordinates occupy one hidden
    column, -1 when they occupy exactly two columns, and zero otherwise.
    Pairs containing the singled-out coordinate are excluded.
    """

    singleton_column = singleton % order
    left_column = pair_left % order
    right_column = pair_right % order
    all_same = (left_column == singleton_column) & (
        right_column == singleton_column
    )
    at_most_two = (
        (left_column == right_column)
        | (left_column == singleton_column)
        | (right_column == singleton_column)
    )
    result = np.where(
        all_same,
        order - 1,
        np.where(at_most_two, -1, 0),
    ).astype(np.int16)
    result[
        (pair_left == singleton) | (pair_right == singleton)
    ] = 0
    return result


def orbit_correlation(order: int) -> np.ndarray:
    """Return the exact integer table C[h,x].

    C[h,x] sums U_0(E) U_h(E) over unordered pairs E of XOR x, where
    U=(q-1)w.  Using unordered pairs directly retains the translation
    stabilizer E xor x = E that the earlier reduction lost.
    """

    dimension = order * order
    if order < 2 or order & (order - 1):
        raise ValueError(("order must be a power of two", order))
    pair_left, pair_right = pair_arrays(dimension)
    pair_xor = pair_left ^ pair_right
    base = scaled_endpoint_weights(
        order, 0, pair_left, pair_right
    ).astype(np.int32)
    correlation = np.empty((dimension, dimension), dtype=np.int64)

    for displacement in range(dimension):
        displaced = scaled_endpoint_weights(
            order, displacement, pair_left, pair_right
        ).astype(np.int32)
        products = base * displaced
        # Every accumulated integer stays below 2^53 through q=32, so the
        # float64 accumulator used by bincount represents it exactly.
        counts = np.bincount(
            pair_xor, weights=products, minlength=dimension
        )
        integer_counts = counts.astype(np.int64)
        if not np.array_equal(counts, integer_counts):
            raise AssertionError(("nonintegral pair count", order))
        correlation[displacement] = integer_counts
    return correlation


def spectrum_certificate(order: int) -> SpectrumCertificate:
    """Compute the exact spectrum numerators and a high-precision coefficient."""

    dimension = order * order
    pair_count = dimension * (dimension - 1) // 2
    scale = (order - 1) ** 4
    correlation = orbit_correlation(order)

    # S(r)=sum_E U_0(E)U_r(E), and B(h,r) is the Walsh transform in pair XOR.
    endpoint_gram = correlation.sum(axis=1, dtype=np.int64)
    pair_symbol = checked_fwht(correlation, axis=1)
    product_bound = int(np.max(np.abs(endpoint_gram))) * int(
        np.max(np.abs(pair_symbol))
    )
    if product_bound > INT64_MAX:
        raise OverflowError(("int64 kernel product", order, product_bound))
    kernel_numerator = endpoint_gram[:, None] * pair_symbol.T

    spectrum_numerator = checked_fwht(
        checked_fwht(kernel_numerator, axis=0),
        axis=1,
    )
    denominator = scale * dimension * dimension
    if np.any(spectrum_numerator <= 0):
        raise AssertionError(
            ("row Gram must be positive definite", order)
        )

    expected_trace = ((dimension + 2) // 2) ** 2
    if int(spectrum_numerator.sum()) != denominator * expected_trace:
        raise AssertionError(
            (
                "block trace",
                order,
                int(spectrum_numerator.sum()),
                denominator * expected_trace,
            )
        )

    numerators, multiplicities = np.unique(
        spectrum_numerator, return_counts=True
    )
    if int(multiplicities.sum()) != dimension * dimension:
        raise AssertionError(("spectrum multiplicity", order))

    with localcontext() as context:
        context.prec = 50
        decimal_denominator = Decimal(denominator)
        root_sum = sum(
            Decimal(int(multiplicity))
            * (Decimal(int(numerator)) / decimal_denominator).sqrt()
            for numerator, multiplicity in zip(
                numerators, multiplicities, strict=True
            )
        )
        coefficient = root_sum / Decimal(dimension * pair_count)

    return SpectrumCertificate(
        order=order,
        dimension=dimension,
        denominator=denominator,
        numerators=numerators,
        multiplicities=multiplicities,
        coefficient=coefficient,
    )


def main() -> None:
    certificates = [spectrum_certificate(order) for order in (2, 4, 8, 16, 32)]
    expected = {
        2: Decimal("0.47159181589114324"),
        4: Decimal("0.06420087162467479"),
    }
    for certificate in certificates:
        if certificate.order in expected:
            error = abs(
                certificate.coefficient - expected[certificate.order]
            )
            if error > Decimal("3e-16"):
                raise AssertionError(
                    (
                        "direct-Gram comparison",
                        certificate.order,
                        certificate.coefficient,
                    )
                )
        print(
            f"q={certificate.order},N={certificate.dimension},"
            f"classes={len(certificate.numerators)},"
            f"rank={int(certificate.multiplicities.sum())},"
            f"coefficient={certificate.coefficient}"
        )

    q32 = certificates[-1]
    print(
        "q=32 exact spectrum numerator/denominator certificate: "
        f"denominator={q32.denominator}"
    )
    for numerator, multiplicity in zip(
        q32.numerators, q32.multiplicities, strict=True
    ):
        print(f"  numerator={int(numerator)},multiplicity={int(multiplicity)}")


if __name__ == "__main__":
    main()
