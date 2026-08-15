#!/usr/bin/env python3
"""Symmetry-reduced exact-moment q=4 screen for the 354 quarantined entries.

Each complement/reversal orbit is represented once.  For every representative
we draw deterministic physical full supports with nonzero signed-permutation
chain moment, cut each support according to the occurrence split, and evaluate
the complete masked kernel on the resulting row and column configurations.
The kernel entries are exact fractions; only the nuclear-norm optimization is
floating point.  Any optimized value above one is a valid lower-witness lead,
while values below one are screens rather than full-operator upper bounds.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from fractions import Fraction
from itertools import product
from json import dumps
from math import prod
from pathlib import Path
import random

import numpy as np
from scipy.optimize import minimize_scalar

from q64_degree_ten_completion_row_insertion import orbit
from q64_masked_universal_audit import affected_classes
from q64_same_side_whole_link_insertion import subset_disjointness_factor
from signed_permutation_link_moment import chain_moment


ROOT = Path(__file__).resolve().parents[1]
ORDER = 4
DIMENSION = ORDER * ORDER
CLASS_NAMES = (
    "universal_septimic",
    "universal_multicubic",
    "universal_double_cubic",
    "universal_noncubic",
    "recovered_universal",
)

Configuration = tuple[tuple[int, ...], ...]
ProfileSplit = tuple[tuple[int, ...], tuple[int, ...]]


@dataclass(frozen=True)
class TemplateRow:
    class_name: str
    profile: tuple[int, ...]
    split: tuple[int, ...]
    orbit_entries: int
    internally_split_blocks: int
    generic_mask_factor: float
    sampled_rows: int
    sampled_columns: int
    exact_nonzero_entries: int
    exact_distinct_values: int
    exact_maximum_absolute_entry: float
    unmasked_uniform_nuclear: float
    masked_uniform_nuclear: float
    masked_optimized_lower: float
    masked_tangent_upper_on_sample: float
    frank_wolfe_iterations: int


@dataclass(frozen=True)
class TemplateScreen:
    order: int
    dimension: int
    affected_entries: int
    complement_reversal_orbits: int
    class_entries: tuple[int, ...]
    class_orbits: tuple[int, ...]
    samples_requested: int
    templates: tuple[TemplateRow, ...]
    maximum_masked_lower: float
    maximum_sample_tangent_upper: float
    templates_with_lower_above_one: int
    templates_with_sample_upper_above_one: int


def template_inventory() -> tuple[tuple[str, ProfileSplit], ...]:
    """Return one canonical representative per complement/reversal orbit."""

    result = []
    for class_name, entries in zip(
        CLASS_NAMES, affected_classes(), strict=True
    ):
        representatives = sorted({min(orbit(entry)) for entry in entries})
        result.extend((class_name, entry) for entry in representatives)
    if len(result) != 97:
        raise AssertionError(("unexpected affected orbit count", len(result)))
    return tuple(result)


def symmetric_difference(
    left: tuple[int, ...], right: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(sorted(set(left).symmetric_difference(right)))


def completed_entry(row: Configuration, column: Configuration) -> Fraction:
    """Return the completed cross-Gram character kernel entry."""

    supports = tuple(
        symmetric_difference(left, right)
        for left, right in zip(row, column, strict=True)
    )
    return chain_moment(ORDER, supports)


def masked_entry(row: Configuration, column: Configuration) -> Fraction:
    """Return the physical occurrence entry with every block mask imposed."""

    if any(
        set(left).intersection(right)
        for left, right in zip(row, column, strict=True)
    ):
        return Fraction(0)
    return completed_entry(row, column)


def random_nonzero_cut(
    profile: tuple[int, ...],
    split: tuple[int, ...],
    rng: random.Random,
) -> tuple[Configuration, Configuration]:
    """Draw one disjoint physical cut whose full-support moment is nonzero."""

    for _ in range(50_000):
        supports = tuple(
            tuple(sorted(rng.sample(range(DIMENSION), degree)))
            for degree in profile
        )
        if not chain_moment(ORDER, supports):
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
    raise RuntimeError(("failed to sample nonzero physical cut", profile, split))


def sampled_configurations(
    entry: ProfileSplit,
    samples: int,
    seed: int,
) -> tuple[tuple[Configuration, ...], tuple[Configuration, ...]]:
    profile, split = entry
    rng = random.Random(seed)
    rows: dict[Configuration, None] = {}
    columns: dict[Configuration, None] = {}
    attempts = 0
    while min(len(rows), len(columns)) < samples:
        row, column = random_nonzero_cut(profile, split, rng)
        rows[row] = None
        columns[column] = None
        attempts += 1
        if attempts > 20 * samples:
            raise RuntimeError(("insufficient distinct sampled cuts", entry))
    return tuple(rows)[:samples], tuple(columns)[:samples]


def exact_matrix(
    rows: tuple[Configuration, ...],
    columns: tuple[Configuration, ...],
    masked: bool,
) -> tuple[np.ndarray, tuple[Fraction, ...]]:
    evaluator = masked_entry if masked else completed_entry
    exact = tuple(
        evaluator(row, column) for row, column in product(rows, columns)
    )
    matrix = np.array(tuple(map(float, exact)), dtype=float).reshape(
        len(rows), len(columns)
    )
    return matrix, exact


def weighted_nuclear(
    kernel: np.ndarray,
    row_law: np.ndarray,
    column_law: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    weighted = (
        np.sqrt(row_law)[:, None]
        * kernel
        * np.sqrt(column_law)[None, :]
    )
    left, singular_values, right = np.linalg.svd(
        weighted, full_matrices=False
    )
    polar = left @ right
    value = float(singular_values.sum())
    row_gradient = np.sum(polar * weighted, axis=1) / (2 * row_law)
    column_gradient = np.sum(polar * weighted, axis=0) / (
        2 * column_law
    )
    return value, row_gradient, column_gradient


def optimize_diagonal_laws(
    kernel: np.ndarray,
    iterations: int = 80,
) -> tuple[float, float, int]:
    """Frank--Wolfe ascent and its concavity tangent upper on one sample."""

    rows, columns = kernel.shape
    row_law = np.full(rows, 1 / rows)
    column_law = np.full(columns, 1 / columns)
    completed_iterations = 0
    for iteration in range(iterations):
        value, row_gradient, column_gradient = weighted_nuclear(
            kernel, row_law, column_law
        )
        tangent_upper = float(
            row_gradient.max() + column_gradient.max()
        )
        if tangent_upper - value <= 1e-9:
            completed_iterations = iteration
            break
        row_vertex = int(np.argmax(row_gradient))
        column_vertex = int(np.argmax(column_gradient))

        def objective(step: float) -> float:
            trial_row = (1 - step) * row_law
            trial_column = (1 - step) * column_law
            trial_row[row_vertex] += step
            trial_column[column_vertex] += step
            return -weighted_nuclear(
                kernel, trial_row, trial_column
            )[0]

        line = minimize_scalar(
            objective,
            bounds=(0.0, 1 - 1e-10),
            method="bounded",
            options={"xatol": 1e-10},
        )
        if not line.success or line.x <= 1e-12:
            completed_iterations = iteration
            break
        row_law *= 1 - line.x
        column_law *= 1 - line.x
        row_law[row_vertex] += line.x
        column_law[column_vertex] += line.x
        completed_iterations = iteration + 1
    value, row_gradient, column_gradient = weighted_nuclear(
        kernel, row_law, column_law
    )
    tangent_upper = float(row_gradient.max() + column_gradient.max())
    return value, tangent_upper, completed_iterations


def entry_mask_factor(entry: ProfileSplit) -> float:
    profile, split = entry
    return prod(
        subset_disjointness_factor(
            min(selected, degree - selected),
            max(selected, degree - selected),
        )
        for degree, selected in zip(profile, split, strict=True)
        if 0 < selected < degree
    )


def screen_template(
    class_name: str,
    entry: ProfileSplit,
    samples: int,
    seed: int,
) -> TemplateRow:
    rows, columns = sampled_configurations(entry, samples, seed)
    unmasked, _ = exact_matrix(rows, columns, masked=False)
    masked, exact = exact_matrix(rows, columns, masked=True)
    row_law = np.full(len(rows), 1 / len(rows))
    column_law = np.full(len(columns), 1 / len(columns))
    unmasked_uniform = weighted_nuclear(
        unmasked, row_law, column_law
    )[0]
    masked_uniform = weighted_nuclear(masked, row_law, column_law)[0]
    optimized, tangent, iterations = optimize_diagonal_laws(masked)
    profile, split = entry
    nonzero = tuple(value for value in exact if value)
    return TemplateRow(
        class_name=class_name,
        profile=profile,
        split=split,
        orbit_entries=len(orbit(entry)),
        internally_split_blocks=sum(
            selected not in (0, degree)
            for degree, selected in zip(profile, split, strict=True)
        ),
        generic_mask_factor=entry_mask_factor(entry),
        sampled_rows=len(rows),
        sampled_columns=len(columns),
        exact_nonzero_entries=len(nonzero),
        exact_distinct_values=len(set(exact)),
        exact_maximum_absolute_entry=float(max(map(abs, exact))),
        unmasked_uniform_nuclear=unmasked_uniform,
        masked_uniform_nuclear=masked_uniform,
        masked_optimized_lower=optimized,
        masked_tangent_upper_on_sample=tangent,
        frank_wolfe_iterations=iterations,
    )


def diagnostic(samples: int = 32) -> TemplateScreen:
    inventory = template_inventory()
    rows = tuple(
        screen_template(
            class_name,
            entry,
            samples,
            seed=10_000 + index,
        )
        for index, (class_name, entry) in enumerate(inventory)
    )
    classes = affected_classes()
    orbit_counts = Counter(name for name, _ in inventory)
    return TemplateScreen(
        order=ORDER,
        dimension=DIMENSION,
        affected_entries=sum(map(len, classes)),
        complement_reversal_orbits=len(inventory),
        class_entries=tuple(map(len, classes)),
        class_orbits=tuple(orbit_counts[name] for name in CLASS_NAMES),
        samples_requested=samples,
        templates=rows,
        maximum_masked_lower=max(row.masked_optimized_lower for row in rows),
        maximum_sample_tangent_upper=max(
            row.masked_tangent_upper_on_sample for row in rows
        ),
        templates_with_lower_above_one=sum(
            row.masked_optimized_lower > 1 for row in rows
        ),
        templates_with_sample_upper_above_one=sum(
            row.masked_tangent_upper_on_sample > 1 for row in rows
        ),
    )


def artifact_text(result: TemplateScreen) -> str:
    payload = {
        "schema": "round4_q4_masked_template_screen_v1",
        "result": asdict(result),
        "evidence_label": (
            "exact q4 signed-permutation moments on deterministic sampled "
            "physical submatrices for all 97 complement/reversal templates; "
            "nuclear norms and concavity tangents are floating diagnostics; "
            "a value below one is not a full-operator theorem"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=32)
    parser.add_argument("--write-artifact", action="store_true")
    arguments = parser.parse_args()
    if arguments.samples < 4:
        raise ValueError("at least four sampled configurations required")
    result = diagnostic(arguments.samples)
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q4_masked_template_screen.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    riskiest = sorted(
        result.templates,
        key=lambda row: row.masked_optimized_lower,
        reverse=True,
    )[:5]
    print(
        "q4 masked template screen: "
        f"entries={result.affected_entries},"
        f"orbits={result.complement_reversal_orbits},"
        f"samples={result.samples_requested},"
        f"maximum_lower={result.maximum_masked_lower:.12g},"
        f"maximum_sample_upper={result.maximum_sample_tangent_upper:.12g},"
        f"lower_above_one={result.templates_with_lower_above_one},"
        f"sample_upper_above_one={result.templates_with_sample_upper_above_one},"
        "riskiest="
        + ";".join(
            f"{row.class_name}:{row.profile}:{row.split}:"
            f"{row.masked_optimized_lower:.8g}"
            for row in riskiest
        )
    )


if __name__ == "__main__":
    main()
