#!/usr/bin/env python3
"""Type-block rank--Frobenius bound for the adjacent mixed-orbit formula.

The exact mixed Fourier reduction has ``N^2`` blocks with rows indexed by
the two selected-pair differences and columns indexed by triple translation
shapes.  The chain symmetry partitions those row and column types into a
finite number of invariant classes.  If ``E_ij`` is the mean squared kernel
energy in a row/column type pair and ``r_ij`` is the smaller type cardinality,
then a candidate Cauchy/rank--Frobenius relaxation assigns

    Phi_ij <= sqrt(r_ij E_ij P_i R_j).

Triangle inequality over type pairs suggests the diagnostic

    Phi <= ||Gamma||_op,    Gamma_ij = sqrt(r_ij E_ij).

The type averages are not proved to control arbitrary within-type physical
laws, so neither the q=4 value nor the sampled moderate-q values are theorem
bounds.  They are route-selection diagnostics for whether this bounded
compound architecture is worth turning into a rigorous model.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import sqrt

import numpy as np

from adjacent_cubic_quintic_mixed_orbit_q4 import (
    build_data,
    combined_link_moment,
)
from adjacent_cubic_quintic_orbit_witness import record_one_link_moment
from opposite_endpoint_orbit_scan import triple_orbit_representatives


@dataclass(frozen=True)
class CompoundTypeBound:
    order: int
    row_types: tuple[tuple, ...]
    column_types: tuple[str, ...]
    row_counts: tuple[int, ...]
    column_counts: tuple[int, ...]
    gamma: np.ndarray
    operator_bound: float
    maximum_energy_spread: float


@dataclass(frozen=True)
class SampledCompoundTypeBound:
    order: int
    samples_per_pair: int
    mean_operator_diagnostic: float
    sampled_max_operator_diagnostic: float
    minimum_energy: float
    maximum_energy: float


def bilinear(left: int, right: int) -> int:
    return int(left & right).bit_count() % 2


def pair_class(left: int, right: int) -> str:
    if left == 0 and right == 0:
        return "00"
    if left == 0:
        return "0*"
    if right == 0:
        return "*0"
    return "**0" if bilinear(left, right) == 0 else "**1"


def row_type(order: int, x: int, y: int) -> tuple[str, str]:
    x_row, x_column = divmod(x, order)
    y_row, y_column = divmod(y, order)
    return pair_class(x_row, y_row), pair_class(x_column, y_column)


def vector_pair_rank(first: int, second: int) -> int:
    if first == 0 and second == 0:
        return 0
    if first == 0 or second == 0 or first == second:
        return 1
    return 2


def vector_pair_kernel(first: int, second: int) -> int | None:
    if vector_pair_rank(first, second) != 1:
        return None
    if first == 0:
        return 1
    if second == 0:
        return 2
    return 3


def column_type(order: int, triple: tuple[int, int, int]) -> str:
    if triple[0] != 0:
        raise ValueError(("canonical triple required", triple))
    first_row, first_column = divmod(triple[1], order)
    second_row, second_column = divmod(triple[2], order)
    row_rank = vector_pair_rank(first_row, second_row)
    column_rank = vector_pair_rank(first_column, second_column)
    if row_rank == column_rank == 1:
        if vector_pair_kernel(first_row, second_row) == vector_pair_kernel(
            first_column,
            second_column,
        ):
            raise AssertionError(("dependent full-cell triple", triple))
        return "11x"
    return f"{row_rank}{column_rank}"


def q4_compound_type_bound() -> CompoundTypeBound:
    data = build_data(4)
    dimension = data.order**2
    row_keys = tuple(
        sorted(
            {
                row_type(data.order, int(x), int(y))
                for x in data.differences
                for y in data.differences
            }
        )
    )
    column_keys = tuple(sorted({column_type(data.order, t) for t in data.triples}))
    row_members = {
        key: [
            (x_index, y_index)
            for x_index, x in enumerate(data.differences)
            for y_index, y in enumerate(data.differences)
            if row_type(data.order, int(x), int(y)) == key
        ]
        for key in row_keys
    }
    column_members = {
        key: [
            index
            for index, triple in enumerate(data.triples)
            if column_type(data.order, triple) == key
        ]
        for key in column_keys
    }
    energies = np.square(data.twisted_spectra).sum(axis=(-2, -1)) / dimension**2
    gamma = np.zeros((len(row_keys), len(column_keys)), dtype=float)
    maximum_spread = 0.0
    for row_index, row_key in enumerate(row_keys):
        for column_index, column_key in enumerate(column_keys):
            values = np.asarray(
                [
                    energies[x_index, y_index, triple_index]
                    for x_index, y_index in row_members[row_key]
                    for triple_index in column_members[column_key]
                ]
            )
            maximum_spread = max(
                maximum_spread,
                float(values.max() - values.min()),
            )
            rank = min(
                len(row_members[row_key]),
                len(column_members[column_key]),
            )
            gamma[row_index, column_index] = sqrt(
                rank * float(values.mean())
            )
    operator = float(np.linalg.svd(gamma, compute_uv=False)[0])
    return CompoundTypeBound(
        order=data.order,
        row_types=row_keys,
        column_types=column_keys,
        row_counts=tuple(len(row_members[key]) for key in row_keys),
        column_counts=tuple(len(column_members[key]) for key in column_keys),
        gamma=gamma,
        operator_bound=operator,
        maximum_energy_spread=maximum_spread,
    )


def kernel_squared_energy(
    order: int,
    x: int,
    y: int,
    triple: tuple[int, int, int],
) -> float:
    """Return sum_{s,t} |M_13 M_35|^2 for one orbit-type atom."""

    dimension = order**2
    total = 0.0
    for s in range(dimension):
        if s in (0, x):
            continue
        cubic = tuple(sorted((0, x, s)))
        middle = record_one_link_moment(order, (0,), cubic)
        if middle == 0:
            continue
        for t in range(dimension):
            shifted = tuple(value ^ t for value in triple)
            if 0 in shifted or y in shifted:
                continue
            quintic = tuple(sorted((0, y) + shifted))
            adjacent = combined_link_moment(order, cubic, quintic)
            total += (middle * adjacent) ** 2
    return total


def sampled_compound_type_bound(
    order: int,
    samples_per_pair: int = 3,
    seed: int = 72031,
) -> SampledCompoundTypeBound:
    """Estimate the type-energy matrix at moderate q.

    The sampled-maximum result is still only a diagnostic: unsampled atoms
    may have larger energy.  It is used to decide whether deriving exact
    general-q type averages is quantitatively promising.
    """

    if order < 4 or order & (order - 1):
        raise ValueError(order)
    dimension = order**2
    differences = tuple(range(1, dimension))
    triples = triple_orbit_representatives(dimension)
    row_members: dict[tuple, list[tuple[int, int]]] = {}
    for x in differences:
        for y in differences:
            row_members.setdefault(row_type(order, x, y), []).append((x, y))
    column_members: dict[str, list[tuple[int, int, int]]] = {}
    for triple in triples:
        column_members.setdefault(column_type(order, triple), []).append(triple)
    row_keys = tuple(sorted(row_members))
    column_keys = tuple(sorted(column_members))
    rng = np.random.default_rng(seed)
    mean_gamma = np.zeros((len(row_keys), len(column_keys)))
    max_gamma = np.zeros_like(mean_gamma)
    minimum_energy = float("inf")
    maximum_energy = 0.0
    for row_index, row_key in enumerate(row_keys):
        for column_index, column_key in enumerate(column_keys):
            rows = row_members[row_key]
            columns = column_members[column_key]
            values = []
            for _ in range(samples_per_pair):
                x, y = rows[int(rng.integers(len(rows)))]
                triple = columns[int(rng.integers(len(columns)))]
                values.append(kernel_squared_energy(order, x, y, triple))
            minimum_energy = min(minimum_energy, min(values))
            maximum_energy = max(maximum_energy, max(values))
            rank = min(len(rows), len(columns))
            mean_gamma[row_index, column_index] = sqrt(
                rank * float(np.mean(values))
            )
            max_gamma[row_index, column_index] = sqrt(rank * max(values))
    return SampledCompoundTypeBound(
        order=order,
        samples_per_pair=samples_per_pair,
        mean_operator_diagnostic=float(
            np.linalg.svd(mean_gamma, compute_uv=False)[0]
        ),
        sampled_max_operator_diagnostic=float(
            np.linalg.svd(max_gamma, compute_uv=False)[0]
        ),
        minimum_energy=minimum_energy,
        maximum_energy=maximum_energy,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-order", type=int)
    parser.add_argument("--samples", type=int, default=3)
    arguments = parser.parse_args()
    if arguments.sample_order:
        sampled = sampled_compound_type_bound(
            arguments.sample_order,
            arguments.samples,
        )
        print(
            "sampled adjacent compound type bound: "
            f"q={sampled.order},samples_per_pair="
            f"{sampled.samples_per_pair},"
            f"mean_operator={sampled.mean_operator_diagnostic:.12g},"
            f"sampled_max_operator="
            f"{sampled.sampled_max_operator_diagnostic:.12g},"
            f"energy_range={sampled.minimum_energy:.12g}:"
            f"{sampled.maximum_energy:.12g}"
        )
        return
    result = q4_compound_type_bound()
    print(
        "adjacent compound type bound: "
        f"q={result.order},row_types={len(result.row_types)},"
        f"column_types={len(result.column_types)},"
        f"row_counts={result.row_counts},"
        f"column_counts={result.column_counts},"
        f"operator_bound={result.operator_bound:.15g},"
        f"maximum_energy_spread={result.maximum_energy_spread:.12g}"
    )
    for row_key, row in zip(result.row_types, result.gamma, strict=True):
        print(
            f"row_type={row_key},gamma="
            + ",".join(f"{value:.12g}" for value in row)
        )


if __name__ == "__main__":
    main()
