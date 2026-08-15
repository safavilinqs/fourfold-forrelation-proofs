#!/usr/bin/env python3
"""Exact obstruction to collectively decoding only two complete active flags.

The accepted active protocol prepares three copies of the folded flag state

    |Psi_x> = (|0>|L_x> + |1>|R_x>)/sqrt(2)

at hard dose two per copy.  This module constructs endpoint ensembles at
N=16 for which even the optimal collective POVM on two complete copies has
Bayes error greater than 1/3.  Tensoring with a fixed F=1 instance lifts the
same obstruction isometrically to N=1024.

The certificate is exact.  It enumerates the N=4 factor ensembles with
F=+/-1/2, forms integer moment matrices, and factors the characteristic
polynomials of sixteen 32-dimensional Gram blocks.  No SDP optimum is used.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from math import prod

import numpy as np
import sympy as sp


BASE_DIMENSION = 4
TENSOR_ENDPOINT_DIMENSION = 16
TARGET_DIMENSION = 1024
MOMENT_DENOMINATOR_PER_COPY = 16
INTEGER_SCALE = 247_808


SquaredSpectrum = tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class TwoFlagCollectiveAudit:
    base_dimension: int
    endpoint_inputs_per_sign: int
    tensor_endpoint_dimension: int
    target_dimension: int
    complete_flags: int
    hard_dose: int
    integer_scale: int
    gram_components: int
    gram_component_size: int
    squared_singular_spectrum: SquaredSpectrum
    nonzero_rank: int
    trace_distance_exact: sp.Expr
    trace_distance: float
    helstrom_error: float
    error_margin_over_one_third: float


def sylvester_sign(dimension: int) -> np.ndarray:
    """Return the unnormalized Sylvester sign matrix."""

    result = np.asarray([[1]], dtype=np.int64)
    while len(result) < dimension:
        result = np.block([[result, result], [result, -result]])
    if result.shape != (dimension, dimension):
        raise ValueError(("Sylvester dimension required", dimension))
    return result


def chain_numerator(
    blocks: np.ndarray,
    hadamard_sign: np.ndarray,
) -> int:
    """Return N^(5/2) F_4 using only integer operations."""

    state = hadamard_sign @ blocks[3]
    state = blocks[2] * state
    state = hadamard_sign @ state
    state = blocks[1] * state
    state = hadamard_sign @ state
    return int(blocks[0] @ state)


@lru_cache(maxsize=1)
def endpoint_factor_data() -> dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Enumerate all N=4 inputs with F=+/-1/2.

    The returned folded-state arrays are scaled by two and hence integral.
    """

    dimension = BASE_DIMENSION
    hadamard_sign = sylvester_sign(dimension)
    folded: dict[int, list[tuple[np.ndarray, np.ndarray, np.ndarray]]] = {
        1: [],
        -1: [],
    }
    for bits in product((-1, 1), repeat=4 * dimension):
        blocks = np.asarray(bits, dtype=np.int64).reshape(4, dimension)
        numerator = chain_numerator(blocks, hadamard_sign)
        if abs(numerator) != 16:
            continue
        sign = 1 if numerator > 0 else -1

        # H_4/2 and |u>/2 imply that 2L and 2R are integral.
        left_numerator = hadamard_sign @ blocks[0]
        right_inner_numerator = hadamard_sign @ blocks[3]
        if np.any(left_numerator % 2):
            raise AssertionError("nonintegral left folded state")
        if np.any(right_inner_numerator % 2):
            raise AssertionError("nonintegral inner right folded state")
        left_scaled = blocks[1] * (left_numerator // 2)
        right_outer_numerator = hadamard_sign @ (
            blocks[2] * (right_inner_numerator // 2)
        )
        if np.any(right_outer_numerator % 2):
            raise AssertionError("nonintegral right folded state")
        right_scaled = right_outer_numerator // 2
        folded[sign].append(
            (blocks.copy(), left_scaled.copy(), right_scaled.copy())
        )

    result: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for sign in (1, -1):
        records = folded[sign]
        result[sign] = (
            np.asarray([record[0] for record in records], dtype=np.int8),
            np.asarray([record[1] for record in records], dtype=np.int8),
            np.asarray([record[2] for record in records], dtype=np.int8),
        )
    return result


def two_copy_moment_numerators(
    sign: int,
) -> tuple[tuple[np.ndarray, ...], ...]:
    """Return exact numerators of the four-by-four path-block moments.

    If V_0=L and V_1=R, path word p=(p_1,p_2) has mode vector
    V_{p_1} tensor V_{p_2}.  Because the stored vectors equal 2V, every
    moment block has common denominator 16 times the ensemble size.
    """

    _, left_scaled, right_scaled = endpoint_factor_data()[sign]
    values = (left_scaled.astype(np.int64), right_scaled.astype(np.int64))
    words = tuple(product((0, 1), repeat=2))
    word_vectors = []
    for first, second in words:
        word_vectors.append(
            np.einsum(
                "ni,nj->nij",
                values[first],
                values[second],
            ).reshape(len(left_scaled), -1)
        )
    return tuple(
        tuple(
            word_vectors[row].T @ word_vectors[column]
            for column in range(4)
        )
        for row in range(4)
    )


def two_flag_integer_cross_block() -> np.ndarray:
    """Return the exact integer-scaled even/odd block of rho_+ - rho_-."""

    moments = {
        sign: two_copy_moment_numerators(sign)
        for sign in (1, -1)
    }
    count = len(endpoint_factor_data()[1][0])
    moment_denominator = MOMENT_DENOMINATOR_PER_COPY * count
    delta_denominator = 8 * moment_denominator**2

    def block(row: int, column: int) -> np.ndarray:
        numerator = sum(
            np.kron(
                moments[sign][row][column],
                moments[sign][row][column],
            )
            - np.kron(
                moments[sign][row][column],
                moments[-sign][row][column],
            )
            for sign in (1, -1)
        )
        scaled = numerator * INTEGER_SCALE
        if np.any(scaled % delta_denominator):
            raise AssertionError(
                ("nonintegral exact cross block", row, column)
            )
        return scaled // delta_denominator

    # Path words are 00, 01, 10, 11.  The hypothesis difference connects
    # the even and odd path-parity sectors and vanishes within each sector.
    return np.block(
        [
            [block(0, 1), block(0, 2)],
            [block(3, 1), block(3, 2)],
        ]
    ).astype(np.int64)


def connected_components(matrix: np.ndarray) -> tuple[tuple[int, ...], ...]:
    """Return components of the nonzero pattern of a symmetric matrix."""

    adjacency = matrix != 0
    visited: set[int] = set()
    components = []
    for root in range(len(matrix)):
        if root in visited:
            continue
        visited.add(root)
        stack = [root]
        component = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor_value in np.flatnonzero(adjacency[current]):
                neighbor = int(neighbor_value)
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        components.append(tuple(sorted(component)))
    return tuple(components)


def exact_squared_singular_spectrum() -> tuple[SquaredSpectrum, int, int]:
    """Factor the exact Gram spectrum component by component."""

    cross = two_flag_integer_cross_block()
    gram = cross @ cross.T
    components = connected_components(gram)
    multiplicities: dict[int, int] = defaultdict(int)
    symbol = sp.Symbol("lambda")
    for component in components:
        block = sp.Matrix(gram[np.ix_(component, component)])
        factors = sp.factor_list(block.charpoly(symbol).as_expr())[1]
        for factor, multiplicity in factors:
            polynomial = sp.Poly(factor, symbol)
            if polynomial.degree() != 1:
                raise AssertionError(("nonlinear Gram factor", factor))
            leading, constant = polynomial.all_coeffs()
            root = sp.Rational(-constant, leading)
            if root.q != 1:
                raise AssertionError(("nonintegral Gram eigenvalue", root))
            multiplicities[int(root)] += int(multiplicity)

    spectrum = tuple(sorted(multiplicities.items(), reverse=True))
    return spectrum, len(components), len(components[0])


EXPECTED_SQUARED_SPECTRUM: SquaredSpectrum = (
    (14_992_384, 1),
    (2_478_080, 12),
    (1_115_136, 6),
    (671_744, 18),
    (184_320, 36),
    (147_456, 18),
    (82_944, 9),
    (0, 412),
)


def trace_distance_from_spectrum(spectrum: SquaredSpectrum) -> sp.Expr:
    """Return ||B||_1, equal to half the trace norm of rho_+ - rho_-."""

    return sp.simplify(
        sum(
            multiplicity * sp.sqrt(eigenvalue)
            for eigenvalue, multiplicity in spectrum
            if eigenvalue
        )
        / INTEGER_SCALE
    )


def fixed_unit_factor() -> np.ndarray:
    """Return one N=4 input with F_4=1."""

    return np.asarray(
        (
            (-1, -1, -1, -1),
            (-1, -1, -1, -1),
            (-1, -1, -1, 1),
            (-1, -1, -1, 1),
        ),
        dtype=np.int64,
    )


def tensor_dimension(*dimensions: int) -> int:
    return prod(dimensions)


def two_flag_collective_audit() -> TwoFlagCollectiveAudit:
    """Run the exact endpoint and Gram-spectrum audit."""

    data = endpoint_factor_data()
    counts = {sign: len(data[sign][0]) for sign in (1, -1)}
    if counts[1] != counts[-1]:
        raise AssertionError(("unbalanced endpoint factor ensembles", counts))

    unit = fixed_unit_factor()
    if chain_numerator(unit, sylvester_sign(4)) != 32:
        raise AssertionError("fixed lift factor does not have F_4=1")
    target = tensor_dimension(TENSOR_ENDPOINT_DIMENSION, 4, 4, 4)
    if target != TARGET_DIMENSION:
        raise AssertionError(("target lift dimension", target))

    spectrum, component_count, component_size = (
        exact_squared_singular_spectrum()
    )
    if spectrum != EXPECTED_SQUARED_SPECTRUM:
        raise AssertionError(("exact Gram spectrum changed", spectrum))
    if sum(multiplicity for _, multiplicity in spectrum) != 512:
        raise AssertionError(("Gram dimension", spectrum))

    distance_exact = trace_distance_from_spectrum(spectrum)
    distance = float(distance_exact)
    error = (1 - distance) / 2
    return TwoFlagCollectiveAudit(
        base_dimension=BASE_DIMENSION,
        endpoint_inputs_per_sign=counts[1],
        tensor_endpoint_dimension=TENSOR_ENDPOINT_DIMENSION,
        target_dimension=target,
        complete_flags=2,
        hard_dose=4,
        integer_scale=INTEGER_SCALE,
        gram_components=component_count,
        gram_component_size=component_size,
        squared_singular_spectrum=spectrum,
        nonzero_rank=sum(
            multiplicity
            for eigenvalue, multiplicity in spectrum
            if eigenvalue
        ),
        trace_distance_exact=distance_exact,
        trace_distance=distance,
        helstrom_error=error,
        error_margin_over_one_third=error - 1 / 3,
    )


def main() -> None:
    result = two_flag_collective_audit()
    print(
        "active two-flag collective obstruction: "
        f"N0={result.tensor_endpoint_dimension},"
        f"N={result.target_dimension},"
        f"factor_inputs_per_sign={result.endpoint_inputs_per_sign},"
        f"flags={result.complete_flags},"
        f"dose={result.hard_dose},"
        f"components={result.gram_components}x"
        f"{result.gram_component_size},"
        f"rank={result.nonzero_rank},"
        f"trace_distance={result.trace_distance:.15g},"
        f"helstrom_error={result.helstrom_error:.15g},"
        f"error_margin={result.error_margin_over_one_third:.15g},"
        f"exact={result.trace_distance_exact}"
    )


if __name__ == "__main__":
    main()
