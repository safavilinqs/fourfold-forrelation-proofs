#!/usr/bin/env python3
"""Evaluate the transposed opposite-endpoint vertical witness.

The previously forced split of profile (3,1,1,5) places endpoint pairs on
one side and the cubic singleton/quintic triple on the other.  The most
expensive remaining split orbit reverses the endpoint support sides while
leaving the two middle singletons in place.  This is a tensor reshuffling,
not a matrix transpose, so the old mixed-orbit value does not automatically
apply.

For the law with all endpoint pairs vertical and all quintic triples in one
hidden column, the row Gram admits three exact reductions:

1. it is block diagonal in h=e xor z;
2. hidden-column translations have zero/nonzero Fourier classes; and
3. a twisted common-row translation has zero/nonzero Fourier classes.

Each final positive-semidefinite block has size q(q-1)^2 instead of the
physical matrix size.  Orders 4 and 8 are fast; order 16 is an optional,
larger exact numerical diagnostic.  No q=32 theorem is claimed here.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import product
from math import comb, sqrt

import numpy as np
from scipy.linalg import qr


@dataclass(frozen=True)
class TransposedVerticalWitness:
    order: int
    dimension: int
    symmetry_block_size: int
    coefficient: float


@dataclass(frozen=True)
class NystromVerticalWitness:
    order: int
    dimension: int
    basis_size: int
    block_ranks: tuple[int, ...]
    coefficient_lower: float


@dataclass(frozen=True)
class CertifiedDominantClassWitness:
    order: int
    dimension: int
    basis_size: int
    retained_rank: int
    common_denominator: int
    contraction_row_upper: float
    spectral_error_upper: float
    computed_trace_lower: float
    rational_trace_lower: int
    coefficient_lower: float


def character_table(order: int) -> np.ndarray:
    return np.asarray(
        [
            [
                -1.0 if int(left & right).bit_count() % 2 else 1.0
                for right in range(order)
            ]
            for left in range(order)
        ]
    )


def state_coordinates(
    order: int,
    indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Decode symmetry-block indices into (x, delta, r) states."""

    q = order
    indices = np.asarray(indices, dtype=np.int64)
    cubic_difference = indices // (q * (q - 1)) + 1
    remainder = indices % (q * (q - 1))
    quintic_difference = remainder // q + 1
    relative_translation = remainder % q
    return (
        cubic_difference.astype(np.int16),
        quintic_difference.astype(np.int16),
        relative_translation.astype(np.int16),
    )


def symmetry_submatrix(
    order: int,
    translation_frequency: int,
    cubic_column_zero: bool,
    quintic_column_zero: bool,
    row_indices: np.ndarray,
    column_indices: np.ndarray,
    *,
    column_chunk_size: int = 16,
) -> np.ndarray:
    """Return a rectangular submatrix without forming the dense full block."""

    q = order
    dimension = q * q
    low_weight = -1 / (q - 1)
    cubic_column_eigenvalue = q - 1 if cubic_column_zero else -1
    quintic_column_eigenvalue = q - 1 if quintic_column_zero else -1
    characters = character_table(q)
    row_x, row_delta, row_relative = state_coordinates(q, row_indices)
    row_x = row_x[:, None]
    row_delta = row_delta[:, None]
    row_relative = row_relative[:, None]
    row_endpoint_xor = np.bitwise_xor(row_relative, row_delta)

    off_column_quintic_sum = (
        2 * low_weight * comb(q - 2, 3)
        + (q - 2) * low_weight**2 * comb(q, 3)
    )
    same_column_quintic_lookup = np.asarray(
        [
            (comb(q - union_size, 3) if q - union_size >= 3 else 0)
            + comb(q, 3) / (q - 1)
            for union_size in range(5)
        ],
        dtype=float,
    )
    result = np.empty((len(row_indices), len(column_indices)))

    for start in range(0, len(column_indices), column_chunk_size):
        stop = min(start + column_chunk_size, len(column_indices))
        selected = np.asarray(column_indices[start:stop])
        column_x, column_delta, column_relative = state_coordinates(q, selected)
        column_x = column_x[None, :]
        column_delta = column_delta[None, :]
        column_relative = column_relative[None, :]

        difference_xor = np.bitwise_xor(row_delta, column_delta)
        zero_character_sum = (difference_xor == 0) * q
        cubic_pair_character = 1 + characters[difference_xor, row_x]
        translated_frequency = np.bitwise_xor(
            translation_frequency,
            column_delta,
        )
        matrix = np.zeros((len(row_indices), stop - start))

        for shift in range(q):
            shifted_relative = np.bitwise_xor(shift, column_relative)
            shifted_endpoint = np.bitwise_xor(
                shifted_relative,
                column_delta,
            )
            intersection = (row_relative == shifted_relative).astype(np.int8)
            intersection += row_relative == shifted_endpoint
            intersection += row_endpoint_xor == shifted_relative
            intersection += row_endpoint_xor == shifted_endpoint
            quintic_factor = (
                same_column_quintic_lookup[4 - intersection]
                + quintic_column_eigenvalue * off_column_quintic_sum
            )

            other_cubic_pair_character = characters[
                difference_xor,
                shift,
            ] * (1 + characters[difference_xor, column_x])
            union_character = (
                cubic_pair_character + other_cubic_pair_character
            )
            first_other_endpoint = np.bitwise_xor(shift, column_x)
            union_character -= (
                (shift == 0) * characters[difference_xor, 0]
            )
            union_character -= (
                (first_other_endpoint == 0)
                * characters[difference_xor, 0]
            )
            union_character -= (
                (row_x == shift) * characters[difference_xor, row_x]
            )
            union_character -= (
                (row_x == first_other_endpoint)
                * characters[difference_xor, row_x]
            )

            same_column_cubic = (
                zero_character_sum
                - union_character
                + (q - 1) * low_weight**2 * zero_character_sum
            )
            off_column_cubic = (
                low_weight * (zero_character_sum - cubic_pair_character)
                + low_weight
                * (zero_character_sum - other_cubic_pair_character)
                + (q - 2) * low_weight**2 * zero_character_sum
            )
            cubic_factor = (
                same_column_cubic
                + cubic_column_eigenvalue * off_column_cubic
            )
            matrix += (
                characters[translated_frequency, shift]
                * quintic_factor
                * cubic_factor
                / dimension**2
            )
        result[:, start:stop] = matrix
    return result


def nystrom_basis_indices(order: int) -> np.ndarray:
    """Choose a full-rank candidate basis from the r=0 state section."""

    q = order
    zero_relative = np.asarray(
        [
            (x - 1) * q * (q - 1) + (delta - 1) * q
            for x in range(1, q)
            for delta in range(1, q)
        ],
        dtype=np.int64,
    )
    reference = symmetry_submatrix(
        q,
        1,
        False,
        False,
        zero_relative,
        zero_relative,
        column_chunk_size=64,
    )
    _, diagonal, pivots = qr(
        reference,
        mode="economic",
        pivoting=True,
    )
    target_rank = q * (q - 1) // 2
    scale = max(1.0, float(abs(diagonal[0, 0])))
    if abs(diagonal[target_rank - 1, target_rank - 1]) < 1e-11 * scale:
        raise AssertionError(("rank-deficient Nystrom basis", q))
    return zero_relative[pivots[:target_rank]]


def nystrom_trace_square_root(
    order: int,
    translation_frequency: int,
    cubic_column_zero: bool,
    quintic_column_zero: bool,
    basis_indices: np.ndarray,
) -> tuple[float, int]:
    """Lower-bound a block trace-square-root by a PSD Nystrom compression."""

    q = order
    size = q * (q - 1) ** 2
    columns = symmetry_submatrix(
        q,
        translation_frequency,
        cubic_column_zero,
        quintic_column_zero,
        np.arange(size, dtype=np.int64),
        basis_indices,
    )
    principal = (columns[basis_indices] + columns[basis_indices].T) / 2
    eigenvalues, eigenvectors = np.linalg.eigh(principal)
    tolerance = 2e-12 * max(1.0, float(eigenvalues[-1]))
    if eigenvalues[0] < -tolerance:
        raise AssertionError(("non-PSD Nystrom principal block", eigenvalues[0]))
    positive = eigenvalues > tolerance
    factor = columns @ (
        eigenvectors[:, positive] / np.sqrt(eigenvalues[positive])
    )
    reduced = (factor.T @ factor)
    reduced = (reduced + reduced.T) / 2
    reduced_eigenvalues = np.linalg.eigvalsh(reduced)
    reduced_tolerance = 2e-9 * max(
        1.0,
        float(reduced_eigenvalues[-1]),
    )
    if reduced_eigenvalues[0] < -reduced_tolerance:
        raise AssertionError(
            ("non-PSD Nystrom reduced block", reduced_eigenvalues[0])
        )
    reduced_eigenvalues[reduced_eigenvalues < 0] = 0
    return (
        float(np.sqrt(reduced_eigenvalues).sum()),
        int(np.count_nonzero(positive)),
    )


def nystrom_vertical_witness(order: int) -> NystromVerticalWitness:
    """Return a deterministic lower witness scalable to q=32."""

    if order < 8 or order & (order - 1):
        raise ValueError(("power-of-two order at least eight required", order))
    q = order
    dimension = q * q
    basis = nystrom_basis_indices(q)
    symmetry_nuclear = 0.0
    ranks = []
    for cubic_zero, quintic_zero in product((True, False), repeat=2):
        column_multiplicity = (
            (1 if cubic_zero else q - 1)
            * (1 if quintic_zero else q - 1)
        )
        for frequency, frequency_multiplicity in ((0, 1), (1, q - 1)):
            trace_root, rank = nystrom_trace_square_root(
                q,
                frequency,
                cubic_zero,
                quintic_zero,
                basis,
            )
            ranks.append(rank)
            symmetry_nuclear += (
                column_multiplicity
                * frequency_multiplicity
                * trace_root
            )

    total_nuclear = dimension * symmetry_nuclear
    triple_count = q * comb(q, 3)
    row_count = dimension**2 * triple_count
    oriented_column_count = dimension**3 * (q - 1) ** 2
    return NystromVerticalWitness(
        order=q,
        dimension=dimension,
        basis_size=len(basis),
        block_ranks=tuple(ranks),
        coefficient_lower=float(
            total_nuclear / sqrt(row_count * oriented_column_count)
        ),
    )


def compensated_matmul(
    left: np.ndarray,
    right: np.ndarray,
    *,
    inner_chunk_size: int = 16,
) -> np.ndarray:
    """Multiply with short dot products and compensated block summation."""

    if left.shape[1] != right.shape[0]:
        raise ValueError(("incompatible matrix shapes", left.shape, right.shape))
    result = np.zeros((left.shape[0], right.shape[1]))
    correction = np.zeros_like(result)
    for start in range(0, left.shape[1], inner_chunk_size):
        stop = min(start + inner_chunk_size, left.shape[1])
        increment = left[:, start:stop] @ right[start:stop]
        adjusted = increment - correction
        updated = result + adjusted
        correction = (updated - result) - adjusted
        result = updated
    return result


def exact_integer_column_gram(
    integer_columns: np.ndarray,
    *,
    digit_bits: int = 16,
) -> np.ndarray:
    """Return C^T C exactly as a Python-integer matrix."""

    count = len(integer_columns)
    base = 1 << digit_bits
    maximum = int(np.max(np.abs(integer_columns)))
    digit_count = max(1, (maximum.bit_length() + digit_bits - 1) // digit_bits)
    signs = np.sign(integer_columns)
    absolute = np.abs(integer_columns)
    digits = tuple(
        (
            signs
            * ((absolute >> (digit_bits * index)) & (base - 1))
        ).astype(float)
        for index in range(digit_count)
    )
    if count * (base - 1) ** 2 >= 2**53:
        raise AssertionError(("digit products are not exactly representable", count))

    result = np.zeros(
        (integer_columns.shape[1], integer_columns.shape[1]),
        dtype=object,
    )
    for left_index, left in enumerate(digits):
        for right_index, right in enumerate(digits):
            product = left.T @ right
            if np.max(np.abs(product - np.rint(product))) != 0:
                raise AssertionError("nonintegral digit product")
            result += (
                product.astype(np.int64).astype(object)
                * base ** (left_index + right_index)
            )
    return result


def certified_dominant_class_witness(
    order: int = 32,
) -> CertifiedDominantClassWitness:
    """Certify one Fourier-class lower witness for the leading split.

    The symmetry-block entries have common denominator
    q^4(q-1)^4.  Recovering their integer numerators makes the selected
    column Gram exact.  A deliberately scaled 480-dimensional factor is
    certified contractive by an absolute-row-sum bound.  Short-dot
    compensated products and explicit roundoff majorants then give a
    lower bound on its trace norm.
    """

    if order != 32:
        raise ValueError("the roundoff certificate is calibrated at q=32")
    q = order
    dimension = q * q
    size = q * (q - 1) ** 2
    denominator = q**4 * (q - 1) ** 4
    basis = nystrom_basis_indices(q)
    columns = symmetry_submatrix(
        q,
        1,
        False,
        False,
        np.arange(size, dtype=np.int64),
        basis,
    )
    scaled_columns = columns * denominator
    integer_columns = np.rint(scaled_columns).astype(np.int64)
    numerator_error = float(
        np.max(np.abs(scaled_columns - integer_columns))
    )
    if numerator_error >= 1 / 64:
        raise AssertionError(("integer-numerator recovery", numerator_error))

    integer_principal = integer_columns[basis]
    if not np.array_equal(integer_principal, integer_principal.T):
        raise AssertionError("nonsymmetric exact principal block")
    principal = integer_principal.astype(float) / denominator
    principal_eigenvalues, principal_eigenvectors = np.linalg.eigh(principal)
    retained = principal_eigenvalues > 1e-5
    if np.count_nonzero(retained) != 480:
        raise AssertionError(("retained certificate rank", np.count_nonzero(retained)))
    safe_factor = 0.9999 * (
        principal_eigenvectors[:, retained]
        / np.sqrt(principal_eigenvalues[retained])
    )

    exact_gram_numerator = exact_integer_column_gram(integer_columns)
    squared_denominator = denominator**2
    gram = np.fromiter(
        (
            float(numerator / squared_denominator)
            for numerator in exact_gram_numerator.flat
        ),
        dtype=float,
        count=len(basis) ** 2,
    ).reshape((len(basis), len(basis)))

    unit_roundoff = 2**-53
    multiplication_factor = 32 * unit_roundoff
    absolute_factor = np.abs(safe_factor)
    principal_entry_error = 1.01 * unit_roundoff * np.abs(principal)
    gram_entry_error = 1.01 * unit_roundoff * np.abs(gram)

    principal_factor = compensated_matmul(principal, safe_factor)
    contraction = compensated_matmul(
        safe_factor.T,
        principal_factor,
    )
    contraction = (contraction + contraction.T) / 2
    principal_factor_error = (
        principal_entry_error @ absolute_factor
        + multiplication_factor
        * (np.abs(principal) @ absolute_factor)
    )
    contraction_error = (
        absolute_factor.T @ principal_factor_error
        + multiplication_factor
        * (absolute_factor.T @ np.abs(principal_factor))
    )
    contraction_row_upper = float(
        np.max(
            np.abs(contraction).sum(axis=1)
            + contraction_error.sum(axis=1)
        )
    )
    if contraction_row_upper >= 1:
        raise AssertionError(("noncontractive safe factor", contraction_row_upper))

    gram_factor = compensated_matmul(gram, safe_factor)
    reduced = compensated_matmul(safe_factor.T, gram_factor)
    reduced = (reduced + reduced.T) / 2
    gram_factor_error = (
        gram_entry_error @ absolute_factor
        + multiplication_factor * (np.abs(gram) @ absolute_factor)
    )
    reduced_error = (
        absolute_factor.T @ gram_factor_error
        + multiplication_factor
        * (absolute_factor.T @ np.abs(gram_factor))
    )
    reduced_eigenvalues, eigenvectors = np.linalg.eigh(reduced)
    diagonalized = compensated_matmul(
        eigenvectors.T,
        compensated_matmul(reduced, eigenvectors),
    )
    off_diagonal = diagonalized - np.diag(np.diag(diagonalized))
    eigenvalue_diagonal_error = float(
        np.max(
            np.abs(
                np.diag(diagonalized) - reduced_eigenvalues
            )
        )
    )
    orthogonality_error = np.linalg.norm(
        eigenvectors.T @ eigenvectors - np.eye(len(eigenvectors)),
        ord="fro",
    )
    eigen_allowance = float(
        eigenvalue_diagonal_error
        + np.linalg.norm(off_diagonal, ord="fro")
        + 2
        * orthogonality_error
        * max(1.0, float(np.max(np.abs(reduced_eigenvalues))))
        + 1e-9
    )
    spectral_error = float(
        np.max(reduced_error.sum(axis=1)) + eigen_allowance
    )
    computed_trace_lower = float(
        np.sqrt(
            np.maximum(reduced_eigenvalues - spectral_error, 0)
        ).sum()
    )
    rational_trace_lower = 194
    if computed_trace_lower <= rational_trace_lower:
        raise AssertionError(
            ("insufficient certified trace lower", computed_trace_lower)
        )

    triple_count = q * comb(q, 3)
    row_count = dimension**2 * triple_count
    oriented_column_count = dimension**3 * (q - 1) ** 2
    # This is just the nonzero-frequency, nonzero/nonzero column class.
    class_multiplicity = (q - 1) ** 3
    coefficient_lower = (
        dimension
        * class_multiplicity
        * rational_trace_lower
        / sqrt(row_count * oriented_column_count)
    )
    return CertifiedDominantClassWitness(
        order=q,
        dimension=dimension,
        basis_size=len(basis),
        retained_rank=int(np.count_nonzero(retained)),
        common_denominator=denominator,
        contraction_row_upper=contraction_row_upper,
        spectral_error_upper=spectral_error,
        computed_trace_lower=computed_trace_lower,
        rational_trace_lower=rational_trace_lower,
        coefficient_lower=float(coefficient_lower),
    )


def symmetry_block(
    order: int,
    translation_frequency: int,
    cubic_column_zero: bool,
    quintic_column_zero: bool,
) -> np.ndarray:
    """Return one exact Fourier block of the oriented-pair column Gram."""

    q = order
    dimension = q * q
    low_weight = -1 / (q - 1)
    cubic_column_eigenvalue = q - 1 if cubic_column_zero else -1
    quintic_column_eigenvalue = q - 1 if quintic_column_zero else -1
    characters = character_table(q)

    # A column state is (x, delta, r): x is the cubic-pair difference,
    # delta the quintic-pair difference, and r their relative translation.
    cubic_difference = np.repeat(np.arange(1, q, dtype=np.int16), q * (q - 1))
    quintic_difference = np.tile(
        np.repeat(np.arange(1, q, dtype=np.int16), q),
        q - 1,
    )
    relative_translation = np.tile(
        np.arange(q, dtype=np.int16),
        (q - 1) ** 2,
    )
    size = len(cubic_difference)

    difference_xor = np.bitwise_xor(
        quintic_difference[:, None],
        quintic_difference[None, :],
    )
    zero_character_sum = (difference_xor == 0) * q
    other_cubic_difference = cubic_difference[None, :]
    other_quintic_difference = quintic_difference[None, :]

    cubic_pair_character = 1 + characters[
        difference_xor,
        cubic_difference[:, None],
    ]
    off_column_quintic_sum = (
        2 * low_weight * comb(q - 2, 3)
        + (q - 2) * low_weight**2 * comb(q, 3)
    )
    same_column_quintic_lookup = np.asarray(
        [
            (comb(q - union_size, 3) if q - union_size >= 3 else 0)
            + comb(q, 3) / (q - 1)
            for union_size in range(5)
        ],
        dtype=float,
    )

    matrix = np.zeros((size, size))
    translated_frequency = np.bitwise_xor(
        translation_frequency,
        other_quintic_difference,
    )

    for shift in range(q):
        shifted_relative = np.bitwise_xor(
            shift,
            relative_translation[None, :],
        )
        shifted_endpoint = np.bitwise_xor(
            shifted_relative,
            other_quintic_difference,
        )
        endpoint = relative_translation[:, None]
        endpoint_xor = np.bitwise_xor(
            relative_translation,
            quintic_difference,
        )[:, None]
        intersection = (endpoint == shifted_relative).astype(np.int8)
        intersection += endpoint == shifted_endpoint
        intersection += endpoint_xor == shifted_relative
        intersection += endpoint_xor == shifted_endpoint
        union_size = 4 - intersection
        quintic_factor = (
            same_column_quintic_lookup[union_size]
            + quintic_column_eigenvalue * off_column_quintic_sum
        )

        other_cubic_pair_character = characters[
            difference_xor,
            shift,
        ] * (
            1
            + characters[
                difference_xor,
                other_cubic_difference,
            ]
        )
        union_character = (
            cubic_pair_character + other_cubic_pair_character
        )
        first_other_endpoint = np.bitwise_xor(
            shift,
            other_cubic_difference,
        )
        # Remove duplicate endpoint characters from the union of the two
        # oriented cubic pairs {0,x} and {shift,shift xor x'}.
        union_character -= (shift == 0) * characters[difference_xor, 0]
        union_character -= (
            first_other_endpoint == 0
        ) * characters[difference_xor, 0]
        union_character -= (
            cubic_difference[:, None] == shift
        ) * characters[difference_xor, cubic_difference[:, None]]
        union_character -= (
            cubic_difference[:, None] == first_other_endpoint
        ) * characters[difference_xor, cubic_difference[:, None]]

        same_column_cubic = (
            zero_character_sum
            - union_character
            + (q - 1) * low_weight**2 * zero_character_sum
        )
        off_column_cubic = (
            low_weight
            * (zero_character_sum - cubic_pair_character)
            + low_weight
            * (zero_character_sum - other_cubic_pair_character)
            + (q - 2) * low_weight**2 * zero_character_sum
        )
        cubic_factor = (
            same_column_cubic
            + cubic_column_eigenvalue * off_column_cubic
        )
        translation_phase = characters[
            translated_frequency,
            shift,
        ]
        matrix += (
            translation_phase
            * quintic_factor
            * cubic_factor
            / dimension**2
        )

    return (matrix + matrix.T) / 2


def trace_square_root(matrix: np.ndarray) -> float:
    eigenvalues = np.linalg.eigvalsh(matrix)
    tolerance = 2e-11 * max(1.0, float(eigenvalues[-1]))
    if eigenvalues[0] < -tolerance:
        raise AssertionError(("non-PSD symmetry block", eigenvalues[0]))
    eigenvalues[eigenvalues < tolerance] = 0
    return float(np.sqrt(eigenvalues).sum())


def transposed_vertical_witness(order: int) -> TransposedVerticalWitness:
    if order < 4 or order & (order - 1):
        raise ValueError(("power-of-two order at least four required", order))
    q = order
    dimension = q * q
    symmetry_nuclear = 0.0
    for cubic_zero, quintic_zero in product((True, False), repeat=2):
        column_multiplicity = (
            (1 if cubic_zero else q - 1)
            * (1 if quintic_zero else q - 1)
        )
        for frequency, frequency_multiplicity in ((0, 1), (1, q - 1)):
            symmetry_nuclear += (
                column_multiplicity
                * frequency_multiplicity
                * trace_square_root(
                    symmetry_block(
                        q,
                        frequency,
                        cubic_zero,
                        quintic_zero,
                    )
                )
            )

    # There are N identical h-blocks.  Oriented endpoint pairs duplicate
    # each physical pair four times and also multiply the column count by
    # four, leaving the normalized coefficient unchanged.
    total_nuclear = dimension * symmetry_nuclear
    triple_count = q * comb(q, 3)
    row_count = dimension**2 * triple_count
    oriented_column_count = dimension**3 * (q - 1) ** 2
    coefficient = total_nuclear / sqrt(row_count * oriented_column_count)
    return TransposedVerticalWitness(
        order=q,
        dimension=dimension,
        symmetry_block_size=q * (q - 1) ** 2,
        coefficient=float(coefficient),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--orders",
        type=int,
        nargs="+",
        default=(4, 8),
        help="powers of two to evaluate; order 16 is substantially larger",
    )
    method = parser.add_mutually_exclusive_group()
    method.add_argument(
        "--nystrom",
        action="store_true",
        help="use the scalable PSD lower compression",
    )
    method.add_argument(
        "--certify-dominant-class",
        action="store_true",
        help="run the exact-numerator q=32 dominant-class certificate",
    )
    arguments = parser.parse_args()
    if arguments.certify_dominant_class:
        result = certified_dominant_class_witness()
        print(
            "certified dominant-class vertical witness: "
            f"q={result.order},"
            f"N={result.dimension},"
            f"basis={result.basis_size},"
            f"retained_rank={result.retained_rank},"
            f"contraction_upper={result.contraction_row_upper:.15g},"
            f"spectral_error={result.spectral_error_upper:.15g},"
            f"computed_trace_lower={result.computed_trace_lower:.15g},"
            f"rational_trace_lower={result.rational_trace_lower},"
            f"coefficient_lower={result.coefficient_lower:.15g}"
        )
        return
    for order in arguments.orders:
        if arguments.nystrom:
            compressed = nystrom_vertical_witness(order)
            print(
                "Nystrom transposed opposite-endpoint witness: "
                f"q={compressed.order},"
                f"N={compressed.dimension},"
                f"basis={compressed.basis_size},"
                f"block_ranks={compressed.block_ranks},"
                f"coefficient_lower={compressed.coefficient_lower:.15g}"
            )
        else:
            result = transposed_vertical_witness(order)
            print(
                "transposed opposite-endpoint vertical witness: "
                f"q={result.order},"
                f"N={result.dimension},"
                f"block_size={result.symmetry_block_size},"
                f"coefficient={result.coefficient:.15g}"
            )


if __name__ == "__main__":
    main()
