#!/usr/bin/env python3
"""Exact masked translation-subspace screens at q=4 and selected q=8.

The coordinate-translation group has ``8 log2(q)`` binary generators across
the four blocks.  A small random linear subspace acts simultaneously on row
and column occurrence configurations.  The complete group-lifted kernel on
that orbit is evaluated exactly, including every distinctness mask.  This
captures coherent multi-mask interference while keeping matrices at most
``2^subspace_dimension`` square.

An optimized value above one would be a physical lower witness.  A value
below one only screens the chosen translation law; it is not an arbitrary-law
upper bound.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from json import dumps
from math import log2
from pathlib import Path
import random

import numpy as np

from q4_masked_template_screen import (
    Configuration,
    optimize_diagonal_laws,
    template_inventory,
)
from signed_permutation_link_moment import chain_moment


ROOT = Path(__file__).resolve().parents[1]
ProfileSplit = tuple[tuple[int, ...], tuple[int, ...]]


@dataclass(frozen=True)
class TranslationRow:
    class_name: str
    profile: tuple[int, ...]
    split: tuple[int, ...]
    selected_trial: int
    group_dimension: int
    group_elements: int
    distinct_row_configurations: int
    distinct_column_configurations: int
    exact_nonzero_entries: int
    exact_distinct_values: int
    exact_maximum_absolute_entry: float
    uniform_nuclear_lower: float
    optimized_nuclear_lower: float
    tangent_upper_on_orbit: float


@dataclass(frozen=True)
class TranslationScreen:
    order: int
    dimension: int
    templates_available: int
    templates_screened: int
    subspace_dimension: int
    trials_per_template: int
    rows: tuple[TranslationRow, ...]
    maximum_optimized_lower: float
    maximum_tangent_upper_on_orbit: float
    templates_above_one: int


def random_nonzero_cut(
    order: int,
    entry: ProfileSplit,
    rng: random.Random,
) -> tuple[Configuration, Configuration]:
    profile, split = entry
    dimension = order * order
    for _ in range(100_000):
        supports = tuple(
            tuple(sorted(rng.sample(range(dimension), degree)))
            for degree in profile
        )
        if not chain_moment(order, supports):
            continue
        row = tuple(
            tuple(sorted(rng.sample(support, selected)))
            for support, selected in zip(supports, split, strict=True)
        )
        column = tuple(
            tuple(value for value in support if value not in set(chosen))
            for support, chosen in zip(supports, row, strict=True)
        )
        return row, column
    raise RuntimeError(("failed to sample nonzero cut", order, entry))


def independent_binary_basis(
    total_bits: int,
    dimension: int,
    rng: random.Random,
) -> tuple[int, ...]:
    pivots: dict[int, int] = {}
    while len(pivots) < dimension:
        value = rng.randrange(1, 1 << total_bits)
        reduced = value
        for pivot in sorted(pivots, reverse=True):
            if (reduced >> pivot) & 1:
                reduced ^= pivots[pivot]
        if not reduced:
            continue
        pivot = reduced.bit_length() - 1
        for old_pivot, old_value in tuple(pivots.items()):
            if (old_value >> pivot) & 1:
                pivots[old_pivot] = old_value ^ reduced
        pivots[pivot] = reduced
    return tuple(pivots[pivot] for pivot in sorted(pivots, reverse=True))


def binary_span(basis: tuple[int, ...]) -> tuple[int, ...]:
    result = [0]
    for vector in basis:
        result += [value ^ vector for value in result]
    return tuple(result)


def translate_configuration(
    order: int,
    configuration: Configuration,
    group_value: int,
) -> Configuration:
    coordinate_bits = 2 * int(log2(order))
    axis_mask = order - 1
    result = []
    for block, support in enumerate(configuration):
        shift = (group_value >> (block * coordinate_bits)) & (
            (1 << coordinate_bits) - 1
        )
        row_shift = shift >> int(log2(order))
        column_shift = shift & axis_mask
        result.append(
            tuple(
                sorted(
                    (((coordinate // order) ^ row_shift) * order)
                    + ((coordinate % order) ^ column_shift)
                    for coordinate in support
                )
            )
        )
    return tuple(result)


def masked_entry(
    order: int,
    row: Configuration,
    column: Configuration,
) -> Fraction:
    if any(
        set(left).intersection(right)
        for left, right in zip(row, column, strict=True)
    ):
        return Fraction(0)
    supports = tuple(
        tuple(sorted(left + right))
        for left, right in zip(row, column, strict=True)
    )
    return chain_moment(order, supports)


def orbit_matrix(
    order: int,
    row: Configuration,
    column: Configuration,
    group: tuple[int, ...],
) -> tuple[np.ndarray, tuple[Fraction, ...], int, int]:
    rows = tuple(translate_configuration(order, row, value) for value in group)
    columns = tuple(
        translate_configuration(order, column, value) for value in group
    )
    exact = tuple(
        masked_entry(order, left, right)
        for left in rows
        for right in columns
    )
    matrix = np.array(tuple(map(float, exact)), dtype=float).reshape(
        len(group), len(group)
    )
    return matrix, exact, len(set(rows)), len(set(columns))


def trial_row(
    order: int,
    class_name: str,
    entry: ProfileSplit,
    subspace_dimension: int,
    trial: int,
    seed: int,
) -> TranslationRow:
    rng = random.Random(seed)
    row, column = random_nonzero_cut(order, entry, rng)
    axis_bits = int(log2(order))
    basis = independent_binary_basis(
        total_bits=8 * axis_bits,
        dimension=subspace_dimension,
        rng=rng,
    )
    group = binary_span(basis)
    matrix, exact, distinct_rows, distinct_columns = orbit_matrix(
        order, row, column, group
    )
    law = np.full(len(group), 1 / len(group))
    uniform = float(
        np.linalg.svd(
            law[:, None] ** 0.5 * matrix * law[None, :] ** 0.5,
            compute_uv=False,
        ).sum()
    )
    optimized, tangent, _ = optimize_diagonal_laws(matrix)
    return TranslationRow(
        class_name=class_name,
        profile=entry[0],
        split=entry[1],
        selected_trial=trial,
        group_dimension=subspace_dimension,
        group_elements=len(group),
        distinct_row_configurations=distinct_rows,
        distinct_column_configurations=distinct_columns,
        exact_nonzero_entries=sum(bool(value) for value in exact),
        exact_distinct_values=len(set(exact)),
        exact_maximum_absolute_entry=float(max(map(abs, exact))),
        uniform_nuclear_lower=uniform,
        optimized_nuclear_lower=optimized,
        tangent_upper_on_orbit=tangent,
    )


def screen_entry(
    order: int,
    class_name: str,
    entry: ProfileSplit,
    subspace_dimension: int,
    trials: int,
    seed: int,
) -> TranslationRow:
    candidates = tuple(
        trial_row(
            order,
            class_name,
            entry,
            subspace_dimension,
            trial,
            seed + trial,
        )
        for trial in range(trials)
    )
    return max(candidates, key=lambda row: row.optimized_nuclear_lower)


def diagnostic(
    order: int = 4,
    subspace_dimension: int = 6,
    trials: int = 3,
    limit: int | None = None,
    selected_indices: tuple[int, ...] | None = None,
) -> TranslationScreen:
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
            subspace_dimension,
            trials,
            seed=90_000 + 100 * index + 10_000 * order,
        )
        for index, (class_name, entry) in enumerate(selected)
    )
    return TranslationScreen(
        order=order,
        dimension=order * order,
        templates_available=len(inventory),
        templates_screened=len(rows),
        subspace_dimension=subspace_dimension,
        trials_per_template=trials,
        rows=rows,
        maximum_optimized_lower=max(
            row.optimized_nuclear_lower for row in rows
        ),
        maximum_tangent_upper_on_orbit=max(
            row.tangent_upper_on_orbit for row in rows
        ),
        templates_above_one=sum(
            row.optimized_nuclear_lower > 1 for row in rows
        ),
    )


def artifact_text(result: TranslationScreen) -> str:
    payload = {
        "schema": "round4_masked_translation_subspace_screen_v1",
        "result": asdict(result),
        "evidence_label": (
            "exact signed-permutation moments on complete deterministic "
            "translation-subspace orbit matrices with every physical mask; "
            "optimized nuclear norms are floating lower witnesses and orbit-"
            "restricted tangents, not full arbitrary-law certificates"
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
    parser.add_argument("--subspace-dimension", type=int, default=6)
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--indices")
    parser.add_argument("--write-artifact", action="store_true")
    arguments = parser.parse_args()
    selected_indices = parse_indices(arguments.indices)
    total_bits = 8 * int(log2(arguments.order))
    if not 1 <= arguments.subspace_dimension <= total_bits:
        raise ValueError(("subspace dimension", arguments.subspace_dimension))
    if arguments.trials < 1:
        raise ValueError("at least one trial required")
    if arguments.limit is not None and arguments.limit < 1:
        raise ValueError("limit must be positive")
    if selected_indices is not None and arguments.limit is not None:
        raise ValueError("use either indices or limit")
    result = diagnostic(
        order=arguments.order,
        subspace_dimension=arguments.subspace_dimension,
        trials=arguments.trials,
        limit=arguments.limit,
        selected_indices=selected_indices,
    )
    if arguments.write_artifact:
        suffix = f"q{arguments.order}_masked_translation_subspace_screen.json"
        path = ROOT / "artifacts" / suffix
        path.write_text(artifact_text(result), encoding="utf-8")
    riskiest = sorted(
        result.rows,
        key=lambda row: row.optimized_nuclear_lower,
        reverse=True,
    )[:8]
    print(
        "masked translation-subspace screen: "
        f"q={result.order},"
        f"templates={result.templates_screened}/{result.templates_available},"
        f"group=2^{result.subspace_dimension},"
        f"trials={result.trials_per_template},"
        f"maximum_lower={result.maximum_optimized_lower:.12g},"
        f"maximum_orbit_upper={result.maximum_tangent_upper_on_orbit:.12g},"
        f"above_one={result.templates_above_one},"
        "riskiest="
        + ";".join(
            f"{row.class_name}:{row.profile}:{row.split}:"
            f"t{row.selected_trial}:{row.optimized_nuclear_lower:.8g}"
            for row in riskiest
        )
    )


if __name__ == "__main__":
    main()
