#!/usr/bin/env python3
"""Optimize mixtures of exact masked translation-orbit shapes.

Exact covariance permits simultaneous translation twirling without loss, but
the twirled laws may mix different translation-orbit shapes.  This screen uses
one common translation subgroup, draws several physical row/column shapes for
each structural template, builds every cross-shape block exactly, and globally
optimizes the two orbit-mass simplexes by concavity-aware Frank--Wolfe ascent.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from math import log2
from pathlib import Path
import random

import numpy as np
from scipy.optimize import minimize_scalar

from masked_translation_subspace_screen import (
    binary_span,
    independent_binary_basis,
    masked_entry,
    random_nonzero_cut,
    translate_configuration,
)
from q4_masked_template_screen import template_inventory, weighted_nuclear


ROOT = Path(__file__).resolve().parents[1]
ProfileSplit = tuple[tuple[int, ...], tuple[int, ...]]


@dataclass(frozen=True)
class MixtureRow:
    class_name: str
    profile: tuple[int, ...]
    split: tuple[int, ...]
    orbit_shapes: int
    group_dimension: int
    group_elements: int
    matrix_dimension: int
    exact_nonzero_entries: int
    exact_distinct_values: int
    uniform_shape_lower: float
    optimized_shape_lower: float
    shape_tangent_upper: float
    row_shape_masses: tuple[float, ...]
    column_shape_masses: tuple[float, ...]
    iterations: int


@dataclass(frozen=True)
class MixtureScreen:
    order: int
    dimension: int
    templates_available: int
    templates_screened: int
    orbit_shapes: int
    group_dimension: int
    rows: tuple[MixtureRow, ...]
    maximum_optimized_lower: float
    maximum_shape_tangent_upper: float
    templates_above_one: int


def mixed_orbit_matrix(
    order: int,
    entry: ProfileSplit,
    orbit_shapes: int,
    group_dimension: int,
    seed: int,
) -> tuple[np.ndarray, int, int]:
    rng = random.Random(seed)
    row_seeds = []
    column_seeds = []
    for _ in range(orbit_shapes):
        row, column = random_nonzero_cut(order, entry, rng)
        row_seeds.append(row)
        column_seeds.append(column)
    axis_bits = int(log2(order))
    basis = independent_binary_basis(
        total_bits=8 * axis_bits,
        dimension=group_dimension,
        rng=rng,
    )
    group = binary_span(basis)
    rows = tuple(
        translate_configuration(order, seed_row, value)
        for seed_row in row_seeds
        for value in group
    )
    columns = tuple(
        translate_configuration(order, seed_column, value)
        for seed_column in column_seeds
        for value in group
    )
    exact = tuple(
        masked_entry(order, row, column)
        for row in rows
        for column in columns
    )
    matrix = np.array(tuple(map(float, exact)), dtype=float).reshape(
        len(rows), len(columns)
    )
    return matrix, sum(bool(value) for value in exact), len(set(exact))


def expanded_law(shape_law: np.ndarray, group_elements: int) -> np.ndarray:
    stable = np.maximum(shape_law, 1e-14)
    stable /= stable.sum()
    return np.repeat(stable / group_elements, group_elements)


def mixture_value_and_gradients(
    matrix: np.ndarray,
    row_shape_law: np.ndarray,
    column_shape_law: np.ndarray,
    group_elements: int,
) -> tuple[float, np.ndarray, np.ndarray]:
    row_law = expanded_law(row_shape_law, group_elements)
    column_law = expanded_law(column_shape_law, group_elements)
    value, row_gradient, column_gradient = weighted_nuclear(
        matrix, row_law, column_law
    )
    row_shape_gradient = row_gradient.reshape(
        len(row_shape_law), group_elements
    ).mean(axis=1)
    column_shape_gradient = column_gradient.reshape(
        len(column_shape_law), group_elements
    ).mean(axis=1)
    return value, row_shape_gradient, column_shape_gradient


def optimize_shape_laws(
    matrix: np.ndarray,
    orbit_shapes: int,
    group_elements: int,
    maximum_iterations: int = 80,
) -> tuple[float, float, np.ndarray, np.ndarray, int]:
    row_law = np.full(orbit_shapes, 1 / orbit_shapes)
    column_law = np.full(orbit_shapes, 1 / orbit_shapes)
    completed = 0
    for iteration in range(maximum_iterations):
        value, row_gradient, column_gradient = mixture_value_and_gradients(
            matrix, row_law, column_law, group_elements
        )
        tangent = float(row_gradient.max() + column_gradient.max())
        if tangent - value <= 1e-10:
            completed = iteration
            break
        row_vertex = int(np.argmax(row_gradient))
        column_vertex = int(np.argmax(column_gradient))

        def objective(step: float) -> float:
            trial_row = (1 - step) * row_law
            trial_column = (1 - step) * column_law
            trial_row[row_vertex] += step
            trial_column[column_vertex] += step
            return -mixture_value_and_gradients(
                matrix,
                trial_row,
                trial_column,
                group_elements,
            )[0]

        line = minimize_scalar(
            objective,
            bounds=(0, 1 - 1e-10),
            method="bounded",
            options={"xatol": 1e-10},
        )
        if not line.success or line.x <= 1e-12:
            completed = iteration
            break
        row_law *= 1 - line.x
        column_law *= 1 - line.x
        row_law[row_vertex] += line.x
        column_law[column_vertex] += line.x
        row_law = np.maximum(row_law, 1e-14)
        column_law = np.maximum(column_law, 1e-14)
        row_law /= row_law.sum()
        column_law /= column_law.sum()
        completed = iteration + 1
    value, row_gradient, column_gradient = mixture_value_and_gradients(
        matrix, row_law, column_law, group_elements
    )
    tangent = float(row_gradient.max() + column_gradient.max())
    return value, tangent, row_law, column_law, completed


def screen_entry(
    order: int,
    class_name: str,
    entry: ProfileSplit,
    orbit_shapes: int,
    group_dimension: int,
    seed: int,
) -> MixtureRow:
    matrix, nonzero, distinct = mixed_orbit_matrix(
        order, entry, orbit_shapes, group_dimension, seed
    )
    group_elements = 1 << group_dimension
    uniform_shape_law = np.full(orbit_shapes, 1 / orbit_shapes)
    uniform = mixture_value_and_gradients(
        matrix,
        uniform_shape_law,
        uniform_shape_law,
        group_elements,
    )[0]
    optimized, tangent, row_law, column_law, iterations = optimize_shape_laws(
        matrix, orbit_shapes, group_elements
    )
    return MixtureRow(
        class_name=class_name,
        profile=entry[0],
        split=entry[1],
        orbit_shapes=orbit_shapes,
        group_dimension=group_dimension,
        group_elements=group_elements,
        matrix_dimension=matrix.shape[0],
        exact_nonzero_entries=nonzero,
        exact_distinct_values=distinct,
        uniform_shape_lower=uniform,
        optimized_shape_lower=optimized,
        shape_tangent_upper=tangent,
        row_shape_masses=tuple(map(float, row_law)),
        column_shape_masses=tuple(map(float, column_law)),
        iterations=iterations,
    )


def diagnostic(
    order: int = 4,
    orbit_shapes: int = 3,
    group_dimension: int = 5,
    limit: int | None = None,
    selected_indices: tuple[int, ...] | None = None,
) -> MixtureScreen:
    inventory = template_inventory()
    if selected_indices is not None:
        selected = tuple(inventory[index] for index in selected_indices)
    elif limit is not None:
        selected = inventory[:limit]
    else:
        selected = inventory
    rows = tuple(
        screen_entry(
            order,
            class_name,
            entry,
            orbit_shapes,
            group_dimension,
            seed=130_000 + 100 * index + 10_000 * order,
        )
        for index, (class_name, entry) in enumerate(selected)
    )
    return MixtureScreen(
        order=order,
        dimension=order * order,
        templates_available=len(inventory),
        templates_screened=len(rows),
        orbit_shapes=orbit_shapes,
        group_dimension=group_dimension,
        rows=rows,
        maximum_optimized_lower=max(row.optimized_shape_lower for row in rows),
        maximum_shape_tangent_upper=max(
            row.shape_tangent_upper for row in rows
        ),
        templates_above_one=sum(row.optimized_shape_lower > 1 for row in rows),
    )


def artifact_text(result: MixtureScreen) -> str:
    payload = {
        "schema": "round4_masked_translation_mixture_screen_v1",
        "result": asdict(result),
        "evidence_label": (
            "exact masked cross-shape translation blocks with floating global "
            "optimization over the selected orbit-shape probability simplexes; "
            "the finite shape sample is diagnostic, not the full twirled theorem"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def parse_indices(value: str | None) -> tuple[int, ...] | None:
    if value is None:
        return None
    return tuple(int(item) for item in value.split(",") if item)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, choices=(4, 8), default=4)
    parser.add_argument("--orbit-shapes", type=int, default=3)
    parser.add_argument("--group-dimension", type=int, default=5)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--indices")
    parser.add_argument("--write-artifact", action="store_true")
    arguments = parser.parse_args()
    selected_indices = parse_indices(arguments.indices)
    if arguments.orbit_shapes < 2:
        raise ValueError("at least two orbit shapes required")
    if arguments.group_dimension < 1:
        raise ValueError("group dimension must be positive")
    if arguments.limit is not None and selected_indices is not None:
        raise ValueError("use either limit or indices")
    result = diagnostic(
        order=arguments.order,
        orbit_shapes=arguments.orbit_shapes,
        group_dimension=arguments.group_dimension,
        limit=arguments.limit,
        selected_indices=selected_indices,
    )
    if arguments.write_artifact:
        suffix = f"q{arguments.order}_masked_translation_mixture_screen.json"
        path = ROOT / "artifacts" / suffix
        path.write_text(artifact_text(result), encoding="utf-8")
    riskiest = sorted(
        result.rows,
        key=lambda row: row.optimized_shape_lower,
        reverse=True,
    )[:8]
    print(
        "masked translation-mixture screen: "
        f"q={result.order},"
        f"templates={result.templates_screened}/{result.templates_available},"
        f"shapes={result.orbit_shapes},"
        f"group=2^{result.group_dimension},"
        f"maximum_lower={result.maximum_optimized_lower:.12g},"
        f"maximum_shape_upper={result.maximum_shape_tangent_upper:.12g},"
        f"above_one={result.templates_above_one},"
        "riskiest="
        + ";".join(
            f"{row.class_name}:{row.profile}:{row.split}:"
            f"{row.optimized_shape_lower:.8g}"
            for row in riskiest
        )
    )


if __name__ == "__main__":
    main()
