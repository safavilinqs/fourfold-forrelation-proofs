#!/usr/bin/env python3
"""Exact full-translation Clifford reduction for masked occurrence kernels.

After gauging the exact row and column translation characters, a lifted
physical kernel has the twisted-convolution form

``K[g,h] = c(g, g xor h) f[g xor h]``.

This module constructs the projective irreducible blocks of that convolution,
including matrix-valued symbols for mixtures of translation-orbit shapes.  It
validates the transform against direct regular matrices, screens one complete
q=4 orbit for all 97 quarantined templates, searches 30 pure shapes in the
leading templates, and evaluates finite mixed-shape concavity tangents.

The exact reduction is theorem infrastructure.  Values below one remain
finite screens until all orbit shapes are bounded.
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

from masked_translation_covariance import translation_cocycle
from masked_translation_covariance import (
    chain_translation_sign,
    link_translation_sign,
    translate_support,
)
from masked_translation_cocycle_inventory import basis_shift
from masked_translation_subspace_screen import (
    masked_entry,
    random_nonzero_cut,
    translate_configuration,
)
from q4_masked_template_screen import template_inventory
from signed_permutation_link_moment import link_moment


ROOT = Path(__file__).resolve().parents[1]


def parity(value: int) -> int:
    return value.bit_count() & 1


def evaluate_form(rows: tuple[int, ...], left: int, right: int) -> int:
    value = 0
    while left:
        bit = (left & -left).bit_length() - 1
        value ^= parity(rows[bit] & right)
        left &= left - 1
    return value


def cocycle_rows(order: int, sizes: tuple[int, int, int, int]) -> tuple[int, ...]:
    d = 8 * int(log2(order))
    basis = tuple(basis_shift(order, bit) for bit in range(d))
    rows = []
    for left in basis:
        row = 0
        for j, right in enumerate(basis):
            if translation_cocycle(order, sizes, left, right) < 0:
                row |= 1 << j
        rows.append(row)
    return tuple(rows)


def alternating_rows(c_rows: tuple[int, ...]) -> tuple[int, ...]:
    d = len(c_rows)
    result = []
    for i in range(d):
        row = 0
        for j in range(d):
            if ((c_rows[i] >> j) ^ (c_rows[j] >> i)) & 1:
                row |= 1 << j
        result.append(row)
    return tuple(result)


def symplectic_basis(
    b_rows: tuple[int, ...],
) -> tuple[tuple[tuple[int, int], ...], tuple[int, ...]]:
    remaining = [1 << bit for bit in range(len(b_rows))]

    def form(left: int, right: int) -> int:
        return evaluate_form(b_rows, left, right)

    pairs: list[tuple[int, int]] = []
    while True:
        found = None
        for i, left in enumerate(remaining):
            for j in range(i + 1, len(remaining)):
                if form(left, remaining[j]):
                    found = (i, j)
                    break
            if found is not None:
                break
        if found is None:
            return tuple(pairs), tuple(remaining)
        i, j = found
        left = remaining[i]
        right = remaining[j]
        survivors = [value for k, value in enumerate(remaining) if k not in (i, j)]
        orthogonal = []
        for value in survivors:
            if form(value, right):
                value ^= left
            if form(value, left):
                value ^= right
            orthogonal.append(value)
        pairs.append((left, right))
        remaining = orthogonal


def fwht_axis(values: np.ndarray, axis: int) -> np.ndarray:
    result = np.swapaxes(np.asarray(values, dtype=complex), axis, -1).copy()
    width = result.shape[-1]
    step = 1
    while step < width:
        blocks = result.reshape(*result.shape[:-1], -1, 2 * step)
        left = blocks[..., :step].copy()
        right = blocks[..., step:].copy()
        blocks[..., :step] = left + right
        blocks[..., step:] = left - right
        step *= 2
    return np.swapaxes(result, axis, -1)


@dataclass(frozen=True)
class System:
    dimension_bits: int
    pairs: tuple[tuple[int, int], ...]
    radical: tuple[int, ...]
    ordered_basis: tuple[int, ...]
    cocycle_rows: tuple[int, ...]
    original_from_coordinates: np.ndarray
    base_phase: np.ndarray

    @property
    def irrep_dimension(self) -> int:
        return 1 << len(self.pairs)

    @property
    def characters(self) -> int:
        return 1 << len(self.radical)


@dataclass(frozen=True)
class PureOrbitRow:
    template_index: int
    class_name: str
    profile: tuple[int, ...]
    split: tuple[int, ...]
    normalized_rank: int
    seed: int
    exact_nonzero_symbol_entries: int
    coefficient: float


@dataclass(frozen=True)
class FocusedPureRow:
    template_index: int
    class_name: str
    profile: tuple[int, ...]
    split: tuple[int, ...]
    normalized_rank: int
    trials: int
    maximum_seed: int
    maximum_nonzero_symbol_entries: int
    maximum_coefficient: float


@dataclass(frozen=True)
class MixedOrbitRow:
    template_index: int
    class_name: str
    profile: tuple[int, ...]
    split: tuple[int, ...]
    normalized_rank: int
    seeds: tuple[int, ...]
    exact_nonzero_symbol_entries: int
    row_shape_law: tuple[float, ...]
    column_shape_law: tuple[float, ...]
    coefficient: float
    concavity_tangent_upper: float


@dataclass(frozen=True)
class FullGroupScreen:
    order: int
    group_elements: int
    templates: int
    canonical_pure_rows: tuple[PureOrbitRow, ...]
    canonical_pure_maximum: float
    focused_pure_rows: tuple[FocusedPureRow, ...]
    focused_pure_maximum: float
    mixed_rows: tuple[MixedOrbitRow, ...]
    mixed_maximum_lower: float
    mixed_maximum_tangent_upper: float
    coefficients_above_one: int


def build_system_from_rows(c_rows: tuple[int, ...]) -> System:
    d = len(c_rows)
    pairs, radical = symplectic_basis(alternating_rows(c_rows))
    ordered = tuple(value for pair in pairs for value in pair) + radical
    if len(ordered) != d:
        raise AssertionError(("incomplete basis", len(ordered), d))
    count = 1 << d
    originals = np.zeros(count, dtype=np.uint32)
    phases = np.ones(count, dtype=complex)
    for coordinates in range(1, count):
        bit = coordinates.bit_length() - 1
        bit_value = 1 << bit
        previous = coordinates ^ bit_value
        current = int(originals[previous])
        generator = ordered[bit]
        sign = -1 if evaluate_form(c_rows, current, generator) else 1
        square_root = 1j if evaluate_form(c_rows, generator, generator) else 1
        originals[coordinates] = current ^ generator
        phases[coordinates] = phases[previous] * sign * square_root
    if len(set(map(int, originals))) != count:
        raise AssertionError("symplectic coordinates are not a basis")
    return System(
        dimension_bits=d,
        pairs=pairs,
        radical=radical,
        ordered_basis=ordered,
        cocycle_rows=c_rows,
        original_from_coordinates=originals,
        base_phase=phases,
    )


def transform(system: System, symbol: np.ndarray) -> np.ndarray:
    n = system.irrep_dimension
    chars = system.characters
    arranged = np.empty((n, n, chars), dtype=complex)
    for coordinates, original in enumerate(system.original_from_coordinates):
        x = 0
        y = 0
        for pair in range(len(system.pairs)):
            x |= ((coordinates >> (2 * pair)) & 1) << pair
            y |= ((coordinates >> (2 * pair + 1)) & 1) << pair
        z = coordinates >> (2 * len(system.pairs))
        arranged[x, y, z] = system.base_phase[coordinates] * symbol[int(original)]
    arranged = fwht_axis(arranged, 2)
    arranged = fwht_axis(arranged, 1)
    matrices = np.empty((chars, n, n), dtype=complex)
    for x in range(n):
        for column in range(n):
            matrices[:, column ^ x, column] = arranged[x, column, :]
    return matrices


def explicit_twisted(c_rows: tuple[int, ...], symbol: np.ndarray) -> np.ndarray:
    count = len(symbol)
    result = np.empty((count, count), dtype=complex)
    for g in range(count):
        for h in range(count):
            t = g ^ h
            sign = -1 if evaluate_form(c_rows, g, t) else 1
            result[g, h] = sign * symbol[t]
    return result


def shifts_from_value(order: int, value: int):
    axis_bits = int(log2(order))
    coordinate_bits = 2 * axis_bits
    mask = order - 1
    return tuple(
        (
            (value >> (block * coordinate_bits + axis_bits)) & mask,
            (value >> (block * coordinate_bits)) & mask,
        )
        for block in range(4)
    )


def cocycle_bilinearity_check(trials: int = 400) -> int:
    """Check the basis-table bicharacter against the exact cocycle formula."""

    rng = random.Random(981)
    size_patterns = (
        (0, 0, 0, 0),
        (0, 0, 1, 1),
        (0, 1, 0, 1),
    )
    checks = 0
    for order in (2, 4, 8, 64):
        dimension_bits = 8 * int(log2(order))
        for sizes in size_patterns:
            rows = cocycle_rows(order, sizes)
            for _ in range(trials):
                left = rng.randrange(1 << dimension_bits)
                right = rng.randrange(1 << dimension_bits)
                exact = translation_cocycle(
                    order,
                    sizes,
                    shifts_from_value(order, left),
                    shifts_from_value(order, right),
                )
                predicted = -1 if evaluate_form(rows, left, right) else 1
                if exact != predicted:
                    raise AssertionError(
                        (
                            "translation cocycle is not the basis bicharacter",
                            order,
                            sizes,
                            left,
                            right,
                            exact,
                            predicted,
                        )
                    )
                checks += 1
    return checks


def physical_transform_check() -> tuple[float, int, int]:
    order = 2
    profile = (1, 1, 3, 1)
    split = (0, 1, 1, 0)
    row, column = random_nonzero_cut(order, (profile, split), random.Random(813))
    row_sizes = tuple(map(len, row))
    column_sizes = tuple(map(len, column))
    row_c = cocycle_rows(order, row_sizes)
    column_c = cocycle_rows(order, column_sizes)
    if row_c != column_c:
        raise AssertionError("complementary cocycles differ")
    count = 1 << len(row_c)
    rows = tuple(translate_configuration(order, row, value) for value in range(count))
    columns = tuple(
        translate_configuration(order, column, value) for value in range(count)
    )
    direct = np.array(
        [float(masked_entry(order, left, right)) for left in rows for right in columns]
    ).reshape(count, count)
    row_phase = np.array(
        [
            chain_translation_sign(order, row, shifts_from_value(order, value))
            for value in range(count)
        ]
    )
    column_phase = np.array(
        [
            chain_translation_sign(order, column, shifts_from_value(order, value))
            for value in range(count)
        ]
    )
    base = np.array(
        [float(masked_entry(order, row, columns[value])) for value in range(count)]
    )
    symbol = column_phase * base
    gauged = row_phase[:, None] * direct * column_phase[None, :]
    expected = explicit_twisted(row_c, symbol)
    if not np.array_equal(gauged, expected):
        raise AssertionError(
            ("physical gauge identity", np.max(abs(gauged - expected)))
        )
    system = build_system_from_rows(row_c)
    transformed = transform(system, symbol)
    direct_value = np.linalg.svd(direct / count, compute_uv=False).sum()
    reduced_value = (
        system.irrep_dimension
        * sum(np.linalg.svd(block, compute_uv=False).sum() for block in transformed)
        / count
    )
    if not np.isclose(direct_value, reduced_value, rtol=1e-12, atol=1e-12):
        raise AssertionError(
            ("physical transformed value", direct_value, reduced_value)
        )
    return direct_value, system.irrep_dimension, system.characters


def full_symbol(order: int, row, column) -> tuple[np.ndarray, int]:
    local_count = order * order
    local_shifts = tuple(
        (value // order, value % order) for value in range(local_count)
    )
    completed = []
    for row_support, column_support in zip(row, column, strict=True):
        translated = tuple(
            translate_support(order, column_support, shift) for shift in local_shifts
        )
        completed.append(
            tuple(
                None
                if set(row_support).intersection(support)
                else tuple(sorted(row_support + support))
                for support in translated
            )
        )
    edges = []
    for block in range(3):
        edge = np.zeros((local_count, local_count), dtype=float)
        for left_value, left_support in enumerate(completed[block]):
            if left_support is None:
                continue
            for right_value, right_support in enumerate(completed[block + 1]):
                if right_support is None:
                    continue
                phase = link_translation_sign(
                    order,
                    column[block],
                    column[block + 1],
                    local_shifts[left_value],
                    local_shifts[right_value],
                )
                edge[left_value, right_value] = phase * float(
                    link_moment(order, left_support, right_support)
                )
        edges.append(edge)
    tensor = (
        edges[0][:, :, None, None]
        * edges[1][None, :, :, None]
        * edges[2][None, None, :, :]
    )
    values = tensor.transpose(3, 2, 1, 0).reshape(-1)
    return values, int(np.count_nonzero(values))


def evaluate_pure_orbit(index: int, seed: int) -> PureOrbitRow:
    order = 4
    class_name, entry = template_inventory()[index]
    row, column = random_nonzero_cut(order, entry, random.Random(seed))
    c_rows = cocycle_rows(order, tuple(map(len, row)))
    complement = cocycle_rows(order, tuple(map(len, column)))
    if c_rows != complement:
        raise AssertionError("row and column cocycles differ")
    system = build_system_from_rows(c_rows)
    symbol, nonzero = full_symbol(order, row, column)
    blocks = transform(system, symbol)
    value = (
        system.irrep_dimension
        * sum(np.linalg.svd(block, compute_uv=False).sum() for block in blocks)
        / len(symbol)
    )
    return PureOrbitRow(
        template_index=index,
        class_name=class_name,
        profile=entry[0],
        split=entry[1],
        normalized_rank=2 * len(system.pairs) // int(log2(order)),
        seed=seed,
        exact_nonzero_symbol_entries=nonzero,
        coefficient=float(value),
    )


def mixed_transforms(order: int, entry, shapes: int, seed: int):
    rng = random.Random(seed)
    rows = []
    columns = []
    for _ in range(shapes):
        row, column = random_nonzero_cut(order, entry, rng)
        rows.append(row)
        columns.append(column)
    return mixed_transforms_from_configurations(order, rows, columns)


def mixed_transforms_from_seeds(order: int, entry, seeds):
    pairs = [random_nonzero_cut(order, entry, random.Random(seed)) for seed in seeds]
    return mixed_transforms_from_configurations(
        order,
        [pair[0] for pair in pairs],
        [pair[1] for pair in pairs],
    )


def mixed_transforms_from_configurations(order: int, rows, columns):
    shapes = len(rows)
    if len(columns) != shapes:
        raise ValueError("row and column shape counts differ")
    c_rows = cocycle_rows(order, tuple(map(len, rows[0])))
    system = build_system_from_rows(c_rows)
    blocks = None
    nonzero = 0
    for i, row in enumerate(rows):
        if cocycle_rows(order, tuple(map(len, row))) != c_rows:
            raise AssertionError("row cocycle changed within sector")
        for j, column in enumerate(columns):
            if cocycle_rows(order, tuple(map(len, column))) != c_rows:
                raise AssertionError("complementary cocycle changed within sector")
            symbol, active = full_symbol(order, row, column)
            nonzero += active
            transformed = transform(system, symbol)
            if blocks is None:
                blocks = np.empty((shapes, shapes, *transformed.shape), dtype=complex)
            blocks[i, j] = transformed
    return system, blocks, nonzero


def mixed_value(
    system: System, blocks: np.ndarray, row_law, column_law, gradients=False
):
    shapes = len(row_law)
    n = system.irrep_dimension
    count = 1 << system.dimension_bits
    root = np.sqrt(row_law)[:, None] * np.sqrt(column_law)[None, :]
    weighted = blocks * root[:, :, None, None, None]
    batch = weighted.transpose(2, 0, 3, 1, 4).reshape(
        system.characters, shapes * n, shapes * n
    )
    scale = n / count
    if not gradients:
        singular = np.linalg.svd(batch, compute_uv=False)
        return float(scale * singular.sum())
    left, singular, right = np.linalg.svd(batch, full_matrices=False)
    polar = left @ right
    products = (np.conj(polar) * batch).real
    value = float(scale * singular.sum())
    row_energy = (
        products.sum(axis=2).reshape(system.characters, shapes, n).sum(axis=(0, 2))
    )
    column_energy = (
        products.sum(axis=1).reshape(system.characters, shapes, n).sum(axis=(0, 2))
    )
    row_gradient = scale * row_energy / (2 * row_law)
    column_gradient = scale * column_energy / (2 * column_law)
    return value, row_gradient, column_gradient


def optimize_mixed(system: System, blocks: np.ndarray, iterations: int = 50):
    shapes = blocks.shape[0]
    row = np.full(shapes, 1 / shapes)
    column = np.full(shapes, 1 / shapes)
    for iteration in range(iterations):
        value, row_gradient, column_gradient = mixed_value(
            system, blocks, row, column, gradients=True
        )
        tangent = float(row_gradient.max() + column_gradient.max())
        if tangent - value < 1e-10:
            break
        i = int(np.argmax(row_gradient))
        j = int(np.argmax(column_gradient))

        def objective(step):
            trial_row = (1 - step) * row
            trial_column = (1 - step) * column
            trial_row[i] += step
            trial_column[j] += step
            return -mixed_value(system, blocks, trial_row, trial_column)

        line = minimize_scalar(
            objective,
            bounds=(0, 1 - 1e-10),
            method="bounded",
            options={"xatol": 1e-10},
        )
        if not line.success or line.x < 1e-12:
            break
        row *= 1 - line.x
        column *= 1 - line.x
        row[i] += line.x
        column[j] += line.x
        row = np.maximum(row, 1e-14)
        row /= row.sum()
        column = np.maximum(column, 1e-14)
        column /= column.sum()
    value, row_gradient, column_gradient = mixed_value(
        system, blocks, row, column, gradients=True
    )
    tangent = float(row_gradient.max() + column_gradient.max())
    return value, tangent, row, column, iteration + 1


def evaluate_mixed_orbits(
    index: int,
    seeds: tuple[int, ...],
    row_law: tuple[float, ...],
    column_law: tuple[float, ...],
) -> MixedOrbitRow:
    class_name, entry = template_inventory()[index]
    system, blocks, nonzero = mixed_transforms_from_seeds(4, entry, seeds)
    row = np.asarray(row_law, dtype=float)
    column = np.asarray(column_law, dtype=float)
    if len(row) != len(seeds) or len(column) != len(seeds):
        raise ValueError("mixed law and seed counts differ")
    if row.sum() <= 0 or column.sum() <= 0:
        raise ValueError("mixed shape laws must have unit mass")
    row = row / row.sum()
    column = column / column.sum()
    value, row_gradient, column_gradient = mixed_value(
        system, blocks, row, column, gradients=True
    )
    tangent = float(row_gradient.max() + column_gradient.max())
    return MixedOrbitRow(
        template_index=index,
        class_name=class_name,
        profile=entry[0],
        split=entry[1],
        normalized_rank=2 * len(system.pairs) // 2,
        seeds=seeds,
        exact_nonzero_symbol_entries=nonzero,
        row_shape_law=tuple(map(float, row)),
        column_shape_law=tuple(map(float, column)),
        coefficient=value,
        concavity_tangent_upper=tangent,
    )


def algebra_transform_check() -> int:
    rng = random.Random(124)
    checks = 0
    for d in range(1, 9):
        for _ in range(10):
            rows = tuple(rng.randrange(1 << d) for _ in range(d))
            system = build_system_from_rows(rows)
            symbol = np.array(
                [
                    rng.randrange(-3, 4) + 1j * rng.randrange(-3, 4)
                    for _ in range(1 << d)
                ]
            )
            direct = np.linalg.svd(
                explicit_twisted(rows, symbol), compute_uv=False
            ).sum()
            blocks = transform(system, symbol)
            reduced = system.irrep_dimension * sum(
                np.linalg.svd(block, compute_uv=False).sum() for block in blocks
            )
            if not np.isclose(direct, reduced, rtol=1e-11, atol=1e-9):
                raise AssertionError((d, direct, reduced, system))
            checks += 1
    for d in range(1, 7):
        rows = tuple(rng.randrange(1 << d) for _ in range(d))
        system = build_system_from_rows(rows)
        symbols = np.array(
            [
                [[rng.randrange(-3, 4) for _ in range(1 << d)] for _ in range(2)]
                for _ in range(2)
            ],
            dtype=float,
        )
        blocks = np.array(
            [[transform(system, symbols[i, j]) for j in range(2)] for i in range(2)]
        )
        p = np.array([0.37, 0.63])
        r = np.array([0.58, 0.42])
        direct_blocks = [
            [
                np.sqrt(p[i] * r[j]) * explicit_twisted(rows, symbols[i, j])
                for j in range(2)
            ]
            for i in range(2)
        ]
        direct = np.linalg.svd(np.block(direct_blocks), compute_uv=False).sum() / (
            1 << d
        )
        reduced, gp, gr = mixed_value(system, blocks, p, r, gradients=True)
        if not np.isclose(direct, reduced, rtol=1e-11, atol=1e-10):
            raise AssertionError(("mixed transform", d, direct, reduced))
        epsilon = 1e-6
        trial = p + np.array([epsilon, -epsilon])
        finite = (mixed_value(system, blocks, trial, r) - reduced) / epsilon
        if not np.isclose(finite, gp[0] - gp[1], rtol=2e-4, atol=2e-4):
            raise AssertionError(("mixed row gradient", d, finite, gp))
        trial = r + np.array([epsilon, -epsilon])
        finite = (mixed_value(system, blocks, p, trial) - reduced) / epsilon
        if not np.isclose(finite, gr[0] - gr[1], rtol=2e-4, atol=2e-4):
            raise AssertionError(("mixed column gradient", d, finite, gr))
        checks += 1
    for d in range(1, 7):
        rows = tuple(rng.randrange(1 << d) for _ in range(d))
        system = build_system_from_rows(rows)
        representations = []
        for group_value in range(1 << d):
            delta = np.zeros(1 << d)
            delta[group_value] = 1
            representations.append(transform(system, delta))
        for left in range(1 << d):
            for right in range(1 << d):
                sign = -1 if evaluate_form(rows, left, right) else 1
                error = np.max(
                    np.abs(
                        representations[left] @ representations[right]
                        - sign * representations[left ^ right]
                    )
                )
                if error:
                    raise AssertionError(
                        ("projective multiplication", d, left, right, error)
                    )
        checks += 1
    return checks


FOCUSED_INDICES = (21, 49, 6, 43, 2, 22, 87)
MIXED_CASES = (
    (
        87,
        (117024, 117000, 117001),
        (0.96198689, 0.01023001, 0.02778310),
        (0.89125122, 0.06094480, 0.04780398),
    ),
    (
        22,
        (52024, 52021, 52029, 52015, 52025),
        (0.25880359, 0.55539976, 0.12806589, 0.02774687, 0.02998390),
        (0.34082596, 0.30143085, 0.19654121, 0.04444092, 0.11676106),
    ),
    (
        49,
        (79023, 79010, 79027),
        (0.99999998, 0.00000001, 0.00000001),
        (0.99999998, 0.00000001, 0.00000001),
    ),
)


def diagnostic() -> FullGroupScreen:
    inventory = template_inventory()
    canonical = tuple(
        evaluate_pure_orbit(index, 7000 + index) for index in range(len(inventory))
    )
    focused = []
    for index in FOCUSED_INDICES:
        trials = tuple(
            evaluate_pure_orbit(
                index,
                30000 + 1000 * index + trial,
            )
            for trial in range(30)
        )
        best = max(trials, key=lambda row: row.coefficient)
        focused.append(
            FocusedPureRow(
                template_index=index,
                class_name=best.class_name,
                profile=best.profile,
                split=best.split,
                normalized_rank=best.normalized_rank,
                trials=len(trials),
                maximum_seed=best.seed,
                maximum_nonzero_symbol_entries=(best.exact_nonzero_symbol_entries),
                maximum_coefficient=best.coefficient,
            )
        )
    mixed = tuple(
        evaluate_mixed_orbits(index, seeds, row_law, column_law)
        for index, seeds, row_law, column_law in MIXED_CASES
    )
    canonical_maximum = max(row.coefficient for row in canonical)
    focused_maximum = max(row.maximum_coefficient for row in focused)
    mixed_maximum = max(row.coefficient for row in mixed)
    tangent_maximum = max(row.concavity_tangent_upper for row in mixed)
    return FullGroupScreen(
        order=4,
        group_elements=1 << 16,
        templates=len(canonical),
        canonical_pure_rows=canonical,
        canonical_pure_maximum=canonical_maximum,
        focused_pure_rows=tuple(focused),
        focused_pure_maximum=focused_maximum,
        mixed_rows=mixed,
        mixed_maximum_lower=mixed_maximum,
        mixed_maximum_tangent_upper=tangent_maximum,
        coefficients_above_one=sum(
            value > 1
            for value in (
                *(row.coefficient for row in canonical),
                *(row.maximum_coefficient for row in focused),
                *(row.concavity_tangent_upper for row in mixed),
            )
        ),
    )


def artifact_text(result: FullGroupScreen) -> str:
    payload = {
        "schema": "round4_q4_masked_translation_full_group_screen_v1",
        "result": asdict(result),
        "evidence_label": (
            "exact signed-permutation moments and exact full-translation "
            "projective block reduction; final singular values and "
            "concavity tangents are floating; pure rows and the listed "
            "finite mixed-shape simplexes are screens, not a bound over all "
            "translation-orbit shapes"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-artifact", action="store_true")
    arguments = parser.parse_args()
    algebra_checks = algebra_transform_check()
    bilinear_checks = cocycle_bilinearity_check()
    direct_value, irrep_dimension, characters = physical_transform_check()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q4_masked_translation_full_group_screen.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q4 masked full-translation screen: "
        f"templates={result.templates},"
        f"canonical_max={result.canonical_pure_maximum:.12g},"
        f"focused_max={result.focused_pure_maximum:.12g},"
        f"mixed_lower={result.mixed_maximum_lower:.12g},"
        f"mixed_tangent={result.mixed_maximum_tangent_upper:.12g},"
        f"above_one={result.coefficients_above_one},"
        f"algebra_checks={algebra_checks},"
        f"bilinear_checks={bilinear_checks},"
        f"q2_direct={direct_value:.12g},"
        f"q2_irrep={irrep_dimension},"
        f"q2_characters={characters},"
        "status=full_group_reduction_valid_no_finite_counterexample"
    )
