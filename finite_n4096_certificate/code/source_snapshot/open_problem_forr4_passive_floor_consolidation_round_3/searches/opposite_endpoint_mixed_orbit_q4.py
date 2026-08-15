#!/usr/bin/env python3
"""Exact q=4 mixed-orbit reduction for the opposite endpoint blocker.

For the critical split of ``(3,1,1,5)``, translation twirling reduces the
row law to ``p[x,y]`` on the two selected-pair XORs and the column law to
``r[tau]`` on triple translation shapes.  Pair-translation Fourier
transforms then give the exact objective

    Phi(p,r) = N^-3 sum_{alpha,gamma} ||L_{alpha,gamma}(p,r)||_*,

where rows of L are indexed by ``(x,y)``, columns by ``tau``, and

    L[(x,y),tau] = sqrt(p[x,y] r[tau])
                     a_x(alpha) g_{y,tau}(gamma xor x).

Here ``a`` is the cubic response spectrum and ``g`` is the
representative-independent signed quintic response spectrum.  This script
validates the reduction against pure fixed-orbit blocks and searches for a
mixed-orbit witness.  It is a finite diagnostic, not a q=32 certificate.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from opposite_endpoint_orbit_scan import (
    cubic_response,
    orbit_block,
    quintic_response,
    support_xor,
    triple_orbit_representatives,
    walsh_transform,
)


@dataclass(frozen=True)
class MixedOrbitData:
    order: int
    differences: np.ndarray
    triples: tuple[tuple[int, int, int], ...]
    cubic_spectrum: np.ndarray
    quintic_spectrum: np.ndarray


@dataclass(frozen=True)
class MixedOrbitEvaluation:
    objective: float
    row_gradient: np.ndarray | None
    column_gradient: np.ndarray | None


def character(left: int, right: int) -> int:
    return -1 if int(left & right).bit_count() % 2 else 1


def build_data(order: int = 4, high_only: bool = False) -> MixedOrbitData:
    dimension = order * order
    differences = np.arange(1, dimension, dtype=np.int64)
    triples = triple_orbit_representatives(dimension)
    cubic = np.asarray(
        [
            walsh_transform(
                cubic_response(order, int(difference), high_only)
            )
            for difference in differences
        ]
    )
    quintic = np.empty(
        (len(differences), len(triples), dimension), dtype=float
    )
    for difference_index, difference in enumerate(differences):
        for triple_index, triple in enumerate(triples):
            spectrum = walsh_transform(
                quintic_response(
                    order,
                    int(difference),
                    triple,
                    high_only,
                )
            )
            triple_xor = support_xor(triple)
            quintic[difference_index, triple_index] = spectrum * np.asarray(
                [
                    character(triple_xor, frequency)
                    for frequency in range(dimension)
                ]
            )
    return MixedOrbitData(
        order=order,
        differences=differences,
        triples=triples,
        cubic_spectrum=cubic,
        quintic_spectrum=quintic,
    )


def evaluate(
    data: MixedOrbitData,
    row_law: np.ndarray,
    column_law: np.ndarray,
    gradients: bool = False,
) -> MixedOrbitEvaluation:
    order = data.order
    dimension = order * order
    pair_types = len(data.differences)
    triple_types = len(data.triples)
    row = np.asarray(row_law, dtype=float)
    column = np.asarray(column_law, dtype=float)
    if row.shape != (pair_types, pair_types):
        raise ValueError(("row-law shape", row.shape))
    if column.shape != (triple_types,):
        raise ValueError(("column-law shape", column.shape))
    if np.any(row < 0) or np.any(column < 0):
        raise ValueError("orbit laws must be nonnegative")
    if not np.isclose(row.sum(), 1) or not np.isclose(column.sum(), 1):
        raise ValueError(("orbit-law mass", row.sum(), column.sum()))
    if gradients and (np.any(row <= 0) or np.any(column <= 0)):
        raise ValueError("gradients require strictly positive laws")

    active_rows = np.argwhere(row > 0)
    active_columns = np.flatnonzero(column > 0)
    x_indices = active_rows[:, 0]
    y_indices = active_rows[:, 1]
    x_differences = data.differences[x_indices]
    row_root = np.sqrt(row[x_indices, y_indices])
    column_root = np.sqrt(column[active_columns])
    row_gradient = np.zeros_like(row) if gradients else None
    column_gradient = np.zeros_like(column) if gradients else None
    objective = 0.0

    for alpha in range(dimension):
        cubic = data.cubic_spectrum[x_indices, alpha]
        for gamma in range(dimension):
            frequencies = np.bitwise_xor(gamma, x_differences)
            quintic = data.quintic_spectrum[
                y_indices[:, None],
                active_columns[None, :],
                frequencies[:, None],
            ]
            matrix = (
                (row_root * cubic)[:, None]
                * quintic
                * column_root[None, :]
            )
            if gradients:
                left, singular, right = np.linalg.svd(
                    matrix, full_matrices=False
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

    scale = dimension**-3
    return MixedOrbitEvaluation(
        objective=objective * scale,
        row_gradient=(row_gradient * scale if gradients else None),
        column_gradient=(column_gradient * scale if gradients else None),
    )


def pure_laws(
    data: MixedOrbitData,
    cubic_difference: int,
    quintic_difference: int,
    triple_shape: tuple[int, int, int],
) -> tuple[np.ndarray, np.ndarray]:
    row = np.zeros((len(data.differences), len(data.differences)))
    column = np.zeros(len(data.triples))
    x = int(np.where(data.differences == cubic_difference)[0][0])
    y = int(np.where(data.differences == quintic_difference)[0][0])
    tau = data.triples.index(triple_shape)
    row[x, y] = 1
    column[tau] = 1
    return row, column


def softmax(values: np.ndarray) -> np.ndarray:
    shifted = values - np.max(values)
    weights = np.exp(shifted)
    return weights / weights.sum()


def optimize_one_law(
    data: MixedOrbitData,
    row_law: np.ndarray,
    column_law: np.ndarray,
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
    data: MixedOrbitData, rounds: int, maximum_iterations: int
) -> tuple[np.ndarray, np.ndarray, MixedOrbitEvaluation]:
    pair_types = len(data.differences)
    triple_types = len(data.triples)
    row = np.full((pair_types, pair_types), 1 / pair_types**2)
    column = np.full(triple_types, 1 / triple_types)
    for iteration in range(rounds):
        row, column, row_value = optimize_one_law(
            data, row, column, True, maximum_iterations
        )
        row, column, column_value = optimize_one_law(
            data, row, column, False, maximum_iterations
        )
        print(
            f"alternating round={iteration + 1},"
            f"row_value={row_value:.12g},"
            f"column_value={column_value:.12g},"
            f"row_max={row.max():.12g},"
            f"column_max={column.max():.12g}"
        )
    return row, column, evaluate(data, row, column, gradients=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--optimize", action="store_true")
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--maxiter", type=int, default=20)
    arguments = parser.parse_args()

    data = build_data()
    cases = (
        (4, 4, (0, 1, 4)),
        (4, 4, (0, 4, 5)),
        (12, 12, (0, 3, 12)),
    )
    for cubic_difference, quintic_difference, triple in cases:
        canonical = min(
            tuple(sorted(value ^ shift for value in triple))
            for shift in triple
        )
        row, column = pure_laws(
            data,
            cubic_difference,
            quintic_difference,
            canonical,
        )
        reduced = evaluate(data, row, column).objective
        direct = orbit_block(
            data.order,
            cubic_difference,
            quintic_difference,
            canonical,
        ).normalized_nuclear
        if not np.isclose(reduced, direct, atol=2e-13):
            raise AssertionError(
                ("pure orbit reduction", reduced, direct, triple)
            )
        print(
            "pure mixed-orbit reduction: "
            f"x={cubic_difference},y={quintic_difference},"
            f"triple={canonical},coefficient={reduced:.12g}"
        )

    pair_types = len(data.differences)
    uniform_row = np.full((pair_types, pair_types), 1 / pair_types**2)
    uniform_column = np.full(len(data.triples), 1 / len(data.triples))
    uniform = evaluate(data, uniform_row, uniform_column, gradients=True)
    print(
        "uniform mixed-orbit law: "
        f"coefficient={uniform.objective:.12g},"
        f"row_gradient_max={uniform.row_gradient.max():.12g},"
        f"column_gradient_max={uniform.column_gradient.max():.12g}"
    )

    if arguments.optimize:
        row, column, result = alternating_search(
            data, arguments.rounds, arguments.maxiter
        )
        row_index = np.unravel_index(np.argmax(row), row.shape)
        column_index = int(np.argmax(column))
        print(
            "mixed-orbit witness: "
            f"coefficient={result.objective:.12g},"
            f"row_type=({int(data.differences[row_index[0]])},"
            f"{int(data.differences[row_index[1]])}),"
            f"row_mass={row[row_index]:.12g},"
            f"triple={data.triples[column_index]},"
            f"column_mass={column[column_index]:.12g},"
            f"row_kkt={result.row_gradient.max()-result.objective/2:.6g},"
            f"column_kkt={result.column_gradient.max()-result.objective/2:.6g}"
        )


if __name__ == "__main__":
    main()
