#!/usr/bin/env python3
"""Small-q orbit scan for the opposite cubic--quintic endpoint chain.

This is an exploratory Track A calculation for the
``(3,1,1,5)/(5,1,1,3)`` blocker.  It uses the exact xor-labelled formulas

    M_31(Q,z) = v_3(Q) H_N(xor(Q),z),
    M_51(Q,z) = v_5(Q) H_N(xor(Q),z),

and factors translation-twirled fixed-split blocks into length-N endpoint
responses.  The scan is a lower-bound/structure diagnostic for physical
orbit laws, not an arbitrary-orbit q=32 upper certificate.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import lru_cache, reduce
from itertools import combinations, product
from math import sqrt
from operator import xor

import numpy as np
from scipy import sparse


@dataclass(frozen=True)
class OrbitBlockResult:
    order: int
    cubic_difference: int
    quintic_difference: int
    triple_shape: tuple[int, ...]
    rows: int
    columns: int
    nonzero: int
    rank: int
    normalized_nuclear: float
    normalized_operator: float


@dataclass(frozen=True)
class EndpointOrbitSummary:
    nuclear: float
    operator: float
    rank: int
    nonzero: int


@dataclass(frozen=True)
class FixedOrbitMaximum:
    order: int
    high_only: bool
    cubic_difference: int
    quintic_difference: int
    triple_shape: tuple[int, int, int]
    coefficient: float


@lru_cache(maxsize=None)
def sylvester(dimension: int) -> np.ndarray:
    result = np.asarray([[1.0]])
    while len(result) < dimension:
        result = np.block([[result, result], [result, -result]])
    return result / sqrt(dimension)


def translated_orbit(
    support: tuple[int, ...], dimension: int
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        sorted(
            {
                tuple(sorted(coordinate ^ shift for coordinate in support))
                for shift in range(dimension)
            }
        )
    )


def translated_sequence(
    support: tuple[int, ...], dimension: int
) -> tuple[tuple[int, ...], ...]:
    result = tuple(
        tuple(sorted(coordinate ^ shift for coordinate in support))
        for shift in range(dimension)
    )
    if len(set(result)) != dimension:
        raise ValueError(("support has a translation stabilizer", support))
    return result


def pair_orbit(
    difference: int, dimension: int
) -> tuple[tuple[int, int], ...]:
    if not 0 < difference < dimension:
        raise ValueError(("nonzero pair difference required", difference))
    return tuple(
        sorted(
            {
                tuple(sorted((coordinate, coordinate ^ difference)))
                for coordinate in range(dimension)
            }
        )
    )


def triple_orbit_representatives(
    dimension: int,
) -> tuple[tuple[int, int, int], ...]:
    result = []
    for first, second in combinations(range(1, dimension), 2):
        shape = (0, first, second)
        canonical = min(
            tuple(sorted(value ^ shift for value in shape))
            for shift in shape
        )
        if shape == canonical:
            result.append(shape)
    expected = (dimension - 1) * (dimension - 2) // 6
    if len(result) != expected:
        raise AssertionError(
            ("triple translation orbits", dimension, len(result))
        )
    return tuple(result)


def support_xor(support: tuple[int, ...]) -> int:
    return reduce(xor, support, 0)


def column_counts(
    support: tuple[int, ...], order: int
) -> dict[int, list[int]]:
    result: dict[int, list[int]] = {}
    for coordinate in support:
        row, column = divmod(coordinate, order)
        result.setdefault(column, []).append(row)
    return result


def cubic_weight(support: tuple[int, ...], order: int) -> float:
    if len(support) != 3:
        raise ValueError(support)
    counts = sorted(len(rows) for rows in column_counts(support, order).values())
    if counts == [3]:
        return 1.0
    if counts == [1, 2]:
        return -1 / (order - 1)
    return 0.0


def quintic_weight(support: tuple[int, ...], order: int) -> float:
    """Exact signed-permutation degree-five endpoint amplitude."""

    if len(support) != 5:
        raise ValueError(support)
    columns = column_counts(support, order)
    odd_columns = [
        column for column, rows in columns.items() if len(rows) % 2
    ]
    if len(odd_columns) != 1:
        return 0.0
    even_xors = []
    for column, rows in columns.items():
        if column == odd_columns[0]:
            continue
        even_xors.append(reduce(xor, rows, 0))
    if not even_xors:
        return 1.0
    if len(even_xors) == 1:
        return 1.0 if even_xors[0] == 0 else -1 / (order - 1)
    if len(even_xors) != 2:
        raise AssertionError(("degree-five even columns", support, columns))
    if even_xors[0] == even_xors[1]:
        return -1 / (order - 1)
    return 2 / ((order - 1) * (order - 2))


def endpoint_moment(
    support: tuple[int, ...],
    singleton: int,
    order: int,
    degree: int,
    high_only: bool,
) -> float:
    weight = (
        cubic_weight(support, order)
        if degree == 3
        else quintic_weight(support, order)
    )
    if high_only and not np.isclose(abs(weight), 1.0):
        return 0.0
    return weight * sylvester(order * order)[support_xor(support), singleton]


def endpoint_incidence_matrices(
    order: int,
    cubic_difference: int,
    quintic_difference: int,
    triple_shape: tuple[int, int, int],
    high_only: bool,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the two endpoint matrices in the factored row Gram."""

    dimension = order * order
    cubic_pairs = pair_orbit(cubic_difference, dimension)
    quintic_pairs = pair_orbit(quintic_difference, dimension)
    triples = translated_orbit(triple_shape, dimension)
    cubic = np.zeros((len(cubic_pairs), dimension))
    quintic = np.zeros((len(quintic_pairs), len(triples)))

    for pair_index, pair in enumerate(cubic_pairs):
        pair_set = set(pair)
        for endpoint_singleton in range(dimension):
            if endpoint_singleton in pair_set:
                continue
            support = tuple(sorted(pair + (endpoint_singleton,)))
            weight = cubic_weight(support, order)
            if high_only and not np.isclose(abs(weight), 1.0):
                continue
            cubic[pair_index, endpoint_singleton] = weight

    for pair_index, pair in enumerate(quintic_pairs):
        pair_set = set(pair)
        for triple_index, triple in enumerate(triples):
            if pair_set.intersection(triple):
                continue
            support = tuple(sorted(pair + triple))
            weight = quintic_weight(support, order)
            if high_only and not np.isclose(abs(weight), 1.0):
                continue
            quintic[pair_index, triple_index] = weight
    return cubic, quintic


def trace_square_root_gram(matrix: np.ndarray) -> tuple[float, int, float]:
    gram = matrix @ matrix.T
    eigenvalues = np.linalg.eigvalsh((gram + gram.T) / 2)
    roots = np.sqrt(np.maximum(eigenvalues, 0))
    return (
        float(roots.sum()),
        int(np.count_nonzero(roots > 2e-10)),
        float(roots[-1]),
    )


def walsh_transform(values: np.ndarray) -> np.ndarray:
    result = np.asarray(values, dtype=float).copy()
    width = 1
    while width < len(result):
        blocks = result.reshape(-1, 2 * width)
        left = blocks[:, :width].copy()
        right = blocks[:, width:].copy()
        blocks[:, :width] = left + right
        blocks[:, width:] = left - right
        width *= 2
    return result


def endpoint_orbit_summary(
    response: np.ndarray, difference: int
) -> EndpointOrbitSummary:
    """Nuclear data for a pair-orbit-by-translation incidence matrix."""

    dimension = len(response)
    if not np.allclose(
        response,
        response[np.arange(dimension) ^ difference],
        atol=2e-14,
    ):
        raise AssertionError(("pair stabilizer", difference))
    frequencies = np.arange(dimension)
    allowed = np.asarray(
        [int(difference & int(value)).bit_count() % 2 == 0 for value in frequencies]
    )
    singular = np.abs(walsh_transform(response)[allowed]) / sqrt(2)
    return EndpointOrbitSummary(
        nuclear=float(singular.sum()),
        operator=float(singular.max(initial=0)),
        rank=int(np.count_nonzero(singular > 2e-10)),
        nonzero=(dimension // 2) * int(np.count_nonzero(response)),
    )


def cubic_response(
    order: int, difference: int, high_only: bool
) -> np.ndarray:
    dimension = order * order
    cubic_pair = (0, difference)
    response = np.zeros(dimension)
    for singleton in range(dimension):
        if singleton in cubic_pair:
            continue
        weight = cubic_weight(tuple(sorted(cubic_pair + (singleton,))), order)
        if not high_only or np.isclose(abs(weight), 1.0):
            response[singleton] = weight
    return response


def cubic_response_summary(
    order: int, difference: int, high_only: bool
) -> EndpointOrbitSummary:
    return endpoint_orbit_summary(
        cubic_response(order, difference, high_only), difference
    )


def quintic_response(
    order: int,
    difference: int,
    triple_shape: tuple[int, int, int],
    high_only: bool,
) -> np.ndarray:
    dimension = order * order
    quintic_pair = (0, difference)
    response = np.zeros(dimension)
    for shift, triple in enumerate(
        translated_sequence(triple_shape, dimension)
    ):
        if set(quintic_pair).intersection(triple):
            continue
        weight = quintic_weight(
            tuple(sorted(quintic_pair + triple)), order
        )
        if not high_only or np.isclose(abs(weight), 1.0):
            response[shift] = weight
    return response


def quintic_response_summary(
    order: int,
    difference: int,
    triple_shape: tuple[int, int, int],
    high_only: bool,
) -> EndpointOrbitSummary:
    return endpoint_orbit_summary(
        quintic_response(
            order, difference, triple_shape, high_only
        ),
        difference,
    )


def orbit_block(
    order: int,
    cubic_difference: int,
    quintic_difference: int,
    triple_shape: tuple[int, int, int],
    high_only: bool = False,
) -> OrbitBlockResult:
    """Evaluate a fixed-shape orbit block through its factored row Gram.

    The split has rows ``(A,C,c)`` and columns ``(e,D,b)``.  Since the
    pair XORs of ``A`` and ``C`` and the triple shape of ``D`` are fixed,
    the endpoint Walsh labels do not depend on the selected pair
    translations.  Orthogonality in ``b`` makes the row Gram

        I_N / N^2 tensor (B_3 B_3^T) tensor (B_5 B_5^T).

    Thus its normalized nuclear norm is computed from two endpoint
    incidence matrices with only ``N/2`` rows each.
    """

    dimension = order * order
    cubic = cubic_response_summary(
        order, cubic_difference, high_only
    )
    quintic = quintic_response_summary(
        order, quintic_difference, triple_shape, high_only
    )
    row_count = dimension**3 // 4
    column_count = dimension**3
    normalization = sqrt(row_count * column_count)
    nuclear = cubic.nuclear * quintic.nuclear
    operator = cubic.operator * quintic.operator / dimension
    return OrbitBlockResult(
        order=order,
        cubic_difference=cubic_difference,
        quintic_difference=quintic_difference,
        triple_shape=triple_shape,
        rows=row_count,
        columns=column_count,
        nonzero=(
            cubic.nonzero
            * quintic.nonzero
            * dimension**2
        ),
        rank=cubic.rank * quintic.rank * dimension,
        normalized_nuclear=float(nuclear / normalization),
        normalized_operator=float(operator / normalization),
    )


def aligned_vertical_orbit_coefficient(
    order: int, high_only: bool
) -> float:
    """Closed value for the aligned vertical-pair orbit family."""

    q = order
    leading = q * (q - 2)
    if high_only:
        return 4 * leading**2 / q**6
    return (
        4
        * (leading + 2)
        * (leading + 2 * (q - 2) / (q - 1))
        / q**6
    )


def exhaustive_fixed_orbit_maximum(
    order: int, high_only: bool
) -> FixedOrbitMaximum:
    """Scan all pair differences and triple shapes at a small order."""

    dimension = order * order
    cubic_nuclear, cubic_difference = max(
        (
            cubic_response_summary(order, difference, high_only).nuclear,
            difference,
        )
        for difference in range(1, dimension)
    )
    quintic_nuclear, quintic_difference, triple_shape = max(
        (
            quintic_response_summary(
                order, difference, shape, high_only
            ).nuclear,
            difference,
            shape,
        )
        for difference in range(1, dimension)
        for shape in triple_orbit_representatives(dimension)
    )
    return FixedOrbitMaximum(
        order=order,
        high_only=high_only,
        cubic_difference=cubic_difference,
        quintic_difference=quintic_difference,
        triple_shape=triple_shape,
        coefficient=(
            2 * cubic_nuclear * quintic_nuclear / dimension**3
        ),
    )


def direct_orbit_block(
    order: int,
    cubic_difference: int,
    quintic_difference: int,
    triple_shape: tuple[int, int, int],
    high_only: bool = False,
) -> OrbitBlockResult:
    """Build the split (2,0,1,2) block for two fixed shape orbits."""

    dimension = order * order
    hadamard = sylvester(dimension)
    cubic_pairs = pair_orbit(cubic_difference, dimension)
    quintic_pairs = pair_orbit(quintic_difference, dimension)
    triples = translated_orbit(triple_shape, dimension)

    left_entries = []
    for pair_index, pair in enumerate(cubic_pairs):
        pair_set = set(pair)
        for endpoint_singleton in range(dimension):
            if endpoint_singleton in pair_set:
                continue
            support = tuple(sorted(pair + (endpoint_singleton,)))
            for first_middle in range(dimension):
                value = endpoint_moment(
                    support,
                    first_middle,
                    order,
                    degree=3,
                    high_only=high_only,
                )
                if value:
                    left_entries.append(
                        (
                            pair_index,
                            endpoint_singleton,
                            first_middle,
                            value,
                        )
                    )

    right_entries = []
    for pair_index, pair in enumerate(quintic_pairs):
        pair_set = set(pair)
        for triple_index, triple in enumerate(triples):
            if pair_set.intersection(triple):
                continue
            support = tuple(sorted(pair + triple))
            for second_middle in range(dimension):
                value = endpoint_moment(
                    support,
                    second_middle,
                    order,
                    degree=5,
                    high_only=high_only,
                )
                if value:
                    right_entries.append(
                        (
                            pair_index,
                            triple_index,
                            second_middle,
                            value,
                        )
                    )

    row_count = len(cubic_pairs) * len(quintic_pairs) * dimension
    column_count = dimension * len(triples) * dimension
    rows = []
    columns = []
    data = []
    for (
        cubic_index,
        endpoint_singleton,
        first_middle,
        left_value,
    ), (
        quintic_index,
        triple_index,
        second_middle,
        right_value,
    ) in product(left_entries, right_entries):
        row = (
            (cubic_index * len(quintic_pairs) + quintic_index)
            * dimension
            + second_middle
        )
        column = (
            (endpoint_singleton * len(triples) + triple_index)
            * dimension
            + first_middle
        )
        rows.append(row)
        columns.append(column)
        data.append(
            left_value
            * hadamard[first_middle, second_middle]
            * right_value
        )
    matrix = sparse.coo_matrix(
        (data, (rows, columns)), shape=(row_count, column_count)
    ).tocsr()
    gram = (matrix @ matrix.T).toarray()
    eigenvalues = np.linalg.eigvalsh((gram + gram.T) / 2)
    singular = np.sqrt(np.maximum(eigenvalues, 0))
    normalization = sqrt(row_count * column_count)
    return OrbitBlockResult(
        order=order,
        cubic_difference=cubic_difference,
        quintic_difference=quintic_difference,
        triple_shape=triple_shape,
        rows=row_count,
        columns=column_count,
        nonzero=matrix.nnz,
        rank=int(np.count_nonzero(singular > 2e-10)),
        normalized_nuclear=float(singular.sum() / normalization),
        normalized_operator=float(singular[-1] / normalization),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--exhaustive",
        action="store_true",
        help="also scan every fixed orbit at q=4 and q=8",
    )
    arguments = parser.parse_args()
    # A vertical pair is two cells in one hidden column.  For the quintic
    # high sector, the triple contains a vertical pair of the same xor plus
    # the odd cell in a second column.
    for order in (4, 8, 16, 32):
        cases = (
            (order, order, (0, order, 1)),
            (order, 2 * order, (0, 2 * order, 1)),
            (order, order, (0, order, order + 1)),
        )
        for high_only in (True, False):
            for cubic_difference, quintic_difference, triple in cases:
                result = orbit_block(
                    order,
                    cubic_difference,
                    quintic_difference,
                    triple,
                    high_only=high_only,
                )
                expected = aligned_vertical_orbit_coefficient(
                    order, high_only
                )
                if not np.isclose(
                    result.normalized_nuclear, expected, atol=2e-14
                ):
                    raise AssertionError(
                        ("aligned vertical closed form", result, expected)
                    )
                if order == 4:
                    direct = direct_orbit_block(
                        order,
                        cubic_difference,
                        quintic_difference,
                        triple,
                        high_only=high_only,
                    )
                    exact_fields = (
                        "order",
                        "cubic_difference",
                        "quintic_difference",
                        "triple_shape",
                        "rows",
                        "columns",
                        "nonzero",
                        "rank",
                    )
                    if any(
                        getattr(result, field) != getattr(direct, field)
                        for field in exact_fields
                    ) or not np.allclose(
                        (
                            result.normalized_nuclear,
                            result.normalized_operator,
                        ),
                        (
                            direct.normalized_nuclear,
                            direct.normalized_operator,
                        ),
                        atol=2e-14,
                    ):
                        raise AssertionError(
                            ("factored/direct orbit block", result, direct)
                        )
                print(
                    "opposite endpoint orbit: "
                    f"q={order},"
                    f"sector={'high' if high_only else 'full'},"
                    f"x={cubic_difference},y={quintic_difference},"
                    f"triple={triple},shape={result.rows}x{result.columns},"
                    f"nonzero={result.nonzero},rank={result.rank},"
                    f"normalized_nuclear={result.normalized_nuclear:.12g},"
                    f"normalized_operator={result.normalized_operator:.12g}"
                )
    if arguments.exhaustive:
        for order in (4, 8):
            for high_only in (True, False):
                result = exhaustive_fixed_orbit_maximum(order, high_only)
                print(
                    "exhaustive fixed-orbit maximum: "
                    f"q={order},"
                    f"sector={'high' if high_only else 'full'},"
                    f"x={result.cubic_difference},"
                    f"y={result.quintic_difference},"
                    f"triple={result.triple_shape},"
                    f"coefficient={result.coefficient:.12g}"
                )


if __name__ == "__main__":
    main()
