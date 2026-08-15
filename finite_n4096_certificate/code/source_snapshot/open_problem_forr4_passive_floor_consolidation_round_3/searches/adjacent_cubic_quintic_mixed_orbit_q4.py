#!/usr/bin/env python3
"""Exact twirled mixed-orbit formula for the adjacent cubic--quintic gate.

For the ``(0,1,2,2)`` split of ``(1,1,3,5)``, translation twirling reduces
the row law to ``Q[x,y]`` on the selected cubic- and quintic-pair XORs and
the column law to ``R[tau]`` on complement-triple translation shapes.

After the first singleton link is removed by orthogonality, use oriented
pairs

    C = (u, u xor x),       D = (v, v xor y),

and write the column singleton and triple translation as

    e = u xor s,            T = v xor t xor tau.

If ``k[x,y,tau](s,t)`` is the product of the remaining two link moments,
put

    H(mu,nu) = sum_{s,t} k(s,t) chi(mu,s) chi(nu,t) chi(s,t).

A Walsh transform in ``u,v,e,T`` followed by one Clifford transform splits
the complete weighted matrix into ``N^2`` blocks

    L[p,n][(x,y),tau]
      = sqrt(Q[x,y] R[tau]) H[p xor xor(tau), n] / N^2.

Thus the arbitrary twirled objective is ``sum_{p,n} ||L[p,n]||_*``.  This
script constructs the exact q=4 blocks, validates pure and mixed laws against
direct occurrence matrices, and can search the two probability simplices.
It is a finite diagnostic, not a q=32 contraction theorem.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import sqrt

import numpy as np
from scipy import sparse
from scipy.optimize import minimize

from adjacent_cubic_quintic_orbit_witness import (
    parity_record_size,
    record_one_link_moment,
    record_three_link_moment,
    unnormalized_sylvester,
    xor_values,
)
from opposite_endpoint_orbit_scan import (
    pair_orbit,
    translated_sequence,
    triple_orbit_representatives,
)


@dataclass(frozen=True)
class AdjacentMixedOrbitData:
    order: int
    differences: np.ndarray
    triples: tuple[tuple[int, int, int], ...]
    triple_xors: np.ndarray
    twisted_spectra: np.ndarray


@dataclass(frozen=True)
class AdjacentMixedOrbitEvaluation:
    objective: float
    row_gradient: np.ndarray | None
    column_gradient: np.ndarray | None


def canonical_triple(
    triple: tuple[int, int, int],
) -> tuple[int, int, int]:
    return min(
        tuple(sorted(value ^ shift for value in triple))
        for shift in triple
    )


def combined_link_moment(
    order: int,
    cubic: tuple[int, int, int],
    quintic: tuple[int, ...],
) -> float:
    """Return M_35 on the only two records allowed after M_13."""

    cubic_record = parity_record_size(order, cubic, axis=1)
    quintic_record = parity_record_size(order, quintic, axis=0)
    if cubic_record != quintic_record:
        return 0.0
    if cubic_record == 1:
        return record_one_link_moment(order, cubic, quintic)
    if cubic_record == 3:
        return record_three_link_moment(order, cubic, quintic)
    return 0.0


def build_data(order: int = 4) -> AdjacentMixedOrbitData:
    """Construct every exact q=4 twisted spectrum."""

    if order != 4:
        raise ValueError("the complete mixed-orbit build is calibrated at q=4")
    dimension = order * order
    differences = np.arange(1, dimension, dtype=np.int16)
    triples = triple_orbit_representatives(dimension)
    triple_xors = np.asarray(
        [xor_values(list(triple)) for triple in triples],
        dtype=np.int16,
    )
    walsh = unnormalized_sylvester(dimension).astype(float)
    spectra = np.zeros(
        (
            len(differences),
            len(differences),
            len(triples),
            dimension,
            dimension,
        ),
        dtype=np.float64,
    )
    twist = walsh
    for x_index, x_value in enumerate(differences):
        x = int(x_value)
        middle = np.zeros(dimension, dtype=float)
        cubics: list[tuple[int, int, int] | None] = []
        for s in range(dimension):
            if s in (0, x):
                cubics.append(None)
                continue
            cubic = tuple(sorted((0, x, s)))
            cubics.append(cubic)
            middle[s] = record_one_link_moment(order, (0,), cubic)
        for y_index, y_value in enumerate(differences):
            y = int(y_value)
            for triple_index, triple in enumerate(triples):
                kernel = np.zeros((dimension, dimension), dtype=float)
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
                        kernel[s, t] = middle[s] * adjacent
                spectra[x_index, y_index, triple_index] = (
                    walsh @ (kernel * twist) @ walsh.T
                )
    return AdjacentMixedOrbitData(
        order=order,
        differences=differences,
        triples=triples,
        triple_xors=triple_xors,
        twisted_spectra=spectra,
    )


def evaluate(
    data: AdjacentMixedOrbitData,
    row_law: np.ndarray,
    column_law: np.ndarray,
    *,
    gradients: bool = False,
) -> AdjacentMixedOrbitEvaluation:
    """Evaluate the exact twirled objective and optional simplex gradients."""

    dimension = data.order**2
    row = np.asarray(row_law, dtype=float)
    column = np.asarray(column_law, dtype=float)
    pair_types = len(data.differences)
    if row.shape != (pair_types, pair_types):
        raise ValueError(("row-law shape", row.shape))
    if column.shape != (len(data.triples),):
        raise ValueError(("column-law shape", column.shape))
    if np.any(row < 0) or np.any(column < 0):
        raise ValueError("orbit laws must be nonnegative")
    if not np.isclose(row.sum(), 1) or not np.isclose(column.sum(), 1):
        raise ValueError(("orbit-law mass", row.sum(), column.sum()))
    active_rows = np.argwhere(row > 0)
    active_columns = np.flatnonzero(column > 0)
    x_indices = active_rows[:, 0]
    y_indices = active_rows[:, 1]
    row_root = np.sqrt(row[x_indices, y_indices])
    column_root = np.sqrt(column[active_columns])
    row_gradient = np.zeros_like(row) if gradients else None
    column_gradient = np.zeros_like(column) if gradients else None
    objective = 0.0
    for p in range(dimension):
        mu = np.bitwise_xor(p, data.triple_xors[active_columns])
        for nu in range(dimension):
            core = data.twisted_spectra[
                x_indices[:, None],
                y_indices[:, None],
                active_columns[None, :],
                mu[None, :],
                nu,
            ]
            matrix = row_root[:, None] * core * column_root[None, :]
            if gradients:
                left, singular, right = np.linalg.svd(
                    matrix,
                    full_matrices=False,
                )
                objective += float(singular.sum())
                row_energy = (left * left) @ singular
                column_energy = (right * right).T @ singular
                row_gradient[x_indices, y_indices] += (
                    row_energy / (2 * row[x_indices, y_indices])
                )
                column_gradient[active_columns] += (
                    column_energy / (2 * column[active_columns])
                )
            else:
                objective += float(
                    np.linalg.svd(matrix, compute_uv=False).sum()
                )
    scale = dimension**-2
    return AdjacentMixedOrbitEvaluation(
        objective=objective * scale,
        row_gradient=(row_gradient * scale if gradients else None),
        column_gradient=(column_gradient * scale if gradients else None),
    )


def pure_laws(
    data: AdjacentMixedOrbitData,
    cubic_difference: int,
    quintic_difference: int,
    triple_shape: tuple[int, int, int],
) -> tuple[np.ndarray, np.ndarray]:
    row = np.zeros(
        (len(data.differences), len(data.differences)),
        dtype=float,
    )
    column = np.zeros(len(data.triples), dtype=float)
    x = int(np.where(data.differences == cubic_difference)[0][0])
    y = int(np.where(data.differences == quintic_difference)[0][0])
    tau = data.triples.index(canonical_triple(triple_shape))
    row[x, y] = 1
    column[tau] = 1
    return row, column


def direct_evaluate_small_support(
    data: AdjacentMixedOrbitData,
    row_law: np.ndarray,
    column_law: np.ndarray,
) -> float:
    """Build a direct reduced matrix for a small-support q=4 law."""

    order = data.order
    dimension = order**2
    active_rows = np.argwhere(row_law > 0)
    active_columns = np.flatnonzero(column_law > 0)
    row_blocks = []
    column_blocks = []
    for x_index, y_index in active_rows:
        x = int(data.differences[x_index])
        y = int(data.differences[y_index])
        for cubic_pair in pair_orbit(x, dimension):
            for quintic_pair in pair_orbit(y, dimension):
                row_blocks.append(
                    (
                        int(x_index),
                        int(y_index),
                        cubic_pair,
                        quintic_pair,
                    )
                )
    for triple_index in active_columns:
        for singleton in range(dimension):
            for triple in translated_sequence(
                data.triples[int(triple_index)],
                dimension,
            ):
                column_blocks.append((int(triple_index), singleton, triple))

    row_indices: list[int] = []
    column_indices: list[int] = []
    values: list[float] = []
    row_denominator = (dimension // 2) ** 2
    column_denominator = dimension**2
    for row_index, (x_index, y_index, cubic_pair, quintic_pair) in enumerate(
        row_blocks
    ):
        cubic_set = set(cubic_pair)
        quintic_set = set(quintic_pair)
        row_weight = sqrt(row_law[x_index, y_index] / row_denominator)
        for column_index, (triple_index, singleton, triple) in enumerate(
            column_blocks
        ):
            if singleton in cubic_set or quintic_set.intersection(triple):
                continue
            cubic = tuple(sorted(cubic_pair + (singleton,)))
            middle = record_one_link_moment(order, (0,), cubic)
            if middle == 0:
                continue
            quintic = tuple(sorted(quintic_pair + triple))
            adjacent = combined_link_moment(order, cubic, quintic)
            if adjacent == 0:
                continue
            column_weight = sqrt(
                column_law[triple_index] / column_denominator
            )
            row_indices.append(row_index)
            column_indices.append(column_index)
            values.append(row_weight * middle * adjacent * column_weight)
    matrix = sparse.coo_matrix(
        (values, (row_indices, column_indices)),
        shape=(len(row_blocks), len(column_blocks)),
    ).toarray()
    return float(np.linalg.svd(matrix, compute_uv=False).sum())


def softmax(values: np.ndarray) -> np.ndarray:
    shifted = values - np.max(values)
    weights = np.exp(shifted)
    return weights / weights.sum()


def optimize_one_law(
    data: AdjacentMixedOrbitData,
    row_law: np.ndarray,
    column_law: np.ndarray,
    *,
    optimize_row: bool,
    maximum_iterations: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    initial_law = row_law.ravel() if optimize_row else column_law
    initial = np.log(initial_law)

    def objective(logits: np.ndarray) -> tuple[float, np.ndarray]:
        law = softmax(logits)
        row = law.reshape(row_law.shape) if optimize_row else row_law
        column = column_law if optimize_row else law
        result = evaluate(data, row, column, gradients=True)
        gradient = (
            result.row_gradient.ravel()
            if optimize_row
            else result.column_gradient
        )
        logit_gradient = law * (gradient - float(law @ gradient))
        return -result.objective, -logit_gradient

    result = minimize(
        objective,
        initial,
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": maximum_iterations, "ftol": 1e-13},
    )
    law = softmax(result.x)
    if optimize_row:
        row_law = law.reshape(row_law.shape)
    else:
        column_law = law
    return row_law, column_law, -float(result.fun)


def alternating_search(
    data: AdjacentMixedOrbitData,
    rounds: int,
    maximum_iterations: int,
) -> tuple[np.ndarray, np.ndarray, AdjacentMixedOrbitEvaluation]:
    pair_types = len(data.differences)
    row = np.full((pair_types, pair_types), 1 / pair_types**2)
    column = np.full(len(data.triples), 1 / len(data.triples))
    for _ in range(rounds):
        row, column, _ = optimize_one_law(
            data,
            row,
            column,
            optimize_row=True,
            maximum_iterations=maximum_iterations,
        )
        row, column, _ = optimize_one_law(
            data,
            row,
            column,
            optimize_row=False,
            maximum_iterations=maximum_iterations,
        )
    return row, column, evaluate(data, row, column, gradients=True)


def fixed_column_scan(
    data: AdjacentMixedOrbitData,
    maximum_iterations: int,
) -> tuple[tuple[float, tuple[int, int, int], np.ndarray], ...]:
    """Optimize the row law for each pure triple orbit."""

    pair_types = len(data.differences)
    results = []
    for triple_index, triple in enumerate(data.triples):
        row = np.full((pair_types, pair_types), 1 / pair_types**2)
        column = np.zeros(len(data.triples))
        column[triple_index] = 1
        row, _, objective = optimize_one_law(
            data,
            row,
            column,
            optimize_row=True,
            maximum_iterations=maximum_iterations,
        )
        results.append((objective, triple, row))
    return tuple(sorted(results, key=lambda item: item[0], reverse=True))


def difference_type(order: int, difference: int) -> str:
    row, column = divmod(difference, order)
    if column == 0:
        return "V"
    if row == 0:
        return "H"
    return "D"


def triple_type(
    order: int,
    triple: tuple[int, int, int],
) -> tuple[tuple[int, ...], tuple[int, ...], str]:
    rows = {}
    columns = {}
    for coordinate in triple:
        row, column = divmod(coordinate, order)
        rows[row] = rows.get(row, 0) + 1
        columns[column] = columns.get(column, 0) + 1
    row_partition = tuple(sorted(rows.values(), reverse=True))
    column_partition = tuple(sorted(columns.values(), reverse=True))
    xor_type = difference_type(order, xor_values(list(triple)))
    return row_partition, column_partition, xor_type


def law_summary(
    data: AdjacentMixedOrbitData,
    row_law: np.ndarray,
    column_law: np.ndarray,
) -> tuple[dict[tuple[str, str], float], dict[tuple, float]]:
    row_summary: dict[tuple[str, str], float] = {}
    for x_index, x in enumerate(data.differences):
        for y_index, y in enumerate(data.differences):
            key = (
                difference_type(data.order, int(x)),
                difference_type(data.order, int(y)),
            )
            row_summary[key] = row_summary.get(key, 0.0) + float(
                row_law[x_index, y_index]
            )
    column_summary: dict[tuple, float] = {}
    for index, triple in enumerate(data.triples):
        key = triple_type(data.order, triple)
        column_summary[key] = column_summary.get(key, 0.0) + float(
            column_law[index]
        )
    return row_summary, column_summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search", action="store_true")
    parser.add_argument("--point-scan", action="store_true")
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--iterations", type=int, default=60)
    arguments = parser.parse_args()
    data = build_data()
    row, column = pure_laws(data, 4, 1, (0, 4, 1))
    pure = evaluate(data, row, column).objective
    print(f"adjacent mixed formula pure check: coefficient={pure:.15g}")
    if arguments.point_scan:
        scan = fixed_column_scan(data, arguments.iterations)
        for rank, (objective, triple, optimized_row) in enumerate(scan[:10], 1):
            row_summary, _ = law_summary(
                data,
                optimized_row,
                np.eye(1, len(data.triples), data.triples.index(triple))[0],
            )
            print(
                "fixed-column adjacent search: "
                f"rank={rank},objective={objective:.15g},"
                f"triple={triple},type={triple_type(data.order, triple)},"
                f"row_types={row_summary}"
            )
    if arguments.search:
        row, column, result = alternating_search(
            data,
            arguments.rounds,
            arguments.iterations,
        )
        print(
            "adjacent mixed-orbit search: "
            f"objective={result.objective:.15g},"
            f"row_support={np.count_nonzero(row > 1e-8)},"
            f"column_support={np.count_nonzero(column > 1e-8)},"
            f"row_max={row.max():.12g},"
            f"column_max={column.max():.12g}"
        )
        row_summary, column_summary = law_summary(data, row, column)
        print(
            "row type masses: "
            + ",".join(
                f"{key[0]}{key[1]}={value:.9g}"
                for key, value in sorted(row_summary.items())
            )
        )
        print(
            "column type masses: "
            + ",".join(
                f"{key}={value:.9g}"
                for key, value in sorted(column_summary.items())
            )
        )
        top_rows = np.argsort(row.ravel())[-12:][::-1]
        print(
            "top row atoms: "
            + ",".join(
                f"({int(data.differences[index // len(data.differences)])},"
                f"{int(data.differences[index % len(data.differences)])})="
                f"{row.ravel()[index]:.9g}"
                for index in top_rows
            )
        )
        top_columns = np.argsort(column)[-12:][::-1]
        print(
            "top column atoms: "
            + ",".join(
                f"{data.triples[index]}={column[index]:.9g}"
                for index in top_columns
            )
        )


if __name__ == "__main__":
    main()
