#!/usr/bin/env python3
"""Symmetric adjacent witness selected by the exact q=4 optimizer.

Both selected pair differences are horizontal, so they are nonzero vectors
``x,y`` in the hidden-column group.  Their law is constant on the two
bilinear classes ``<x,y>=0`` and ``<x,y>=1``.  The complement triple is
vertical, and its law is uniform over every translation orbit of three
distinct hidden rows.  The diagonal ``GL(m,2)`` action preserves the
bilinear class and is transitive on the triple shapes, so this is a closed
one-parameter physical family.

The script constructs every exact twisted spectrum at q<=8.  It also derives
a 25-frequency-class closed formula, validates that formula against the
complete q=4 and q=8 constructions, and evaluates the same physical family
through q=32 without storing the full mixed-orbit tensor.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from adjacent_cubic_quintic_mixed_orbit_q4 import combined_link_moment
from adjacent_cubic_quintic_orbit_witness import (
    character,
    parity_record_size,
    record_one_link_moment,
    xor_values,
)
from adjacent_cubic_quintic_structured_witness import twisted_walsh
from opposite_endpoint_orbit_scan import triple_orbit_representatives


@dataclass(frozen=True)
class SymmetricVerticalTripleData:
    order: int
    differences: np.ndarray
    triples: tuple[tuple[int, int, int], ...]
    triple_xors: np.ndarray
    record_one_spectra: np.ndarray
    record_three_spectra: np.ndarray


@dataclass(frozen=True)
class SymmetricVerticalTripleWitness:
    order: int
    orthogonal_mass: float
    orthogonal_pairs: int
    nonorthogonal_pairs: int
    triple_orbits: int
    coefficient: float
    record_one_triangle: float
    record_three_triangle: float


@dataclass(frozen=True)
class ClosedSymmetricVerticalTripleWitness:
    order: int
    orthogonal_mass: float
    orthogonal_pairs: int
    nonorthogonal_pairs: int
    triple_orbits: int
    frequency_classes: int
    coefficient: float


def pair_frequency_class(left: int, right: int) -> str:
    if left == 0 and right == 0:
        return "00"
    if left == 0:
        return "0*"
    if right == 0:
        return "*0"
    return "**0" if bilinear(left, right) == 0 else "**1"


def bilinear(left: int, right: int) -> int:
    return int(left & right).bit_count() % 2


def character_sum_excluding_pair(
    order: int,
    label: int,
    excluded: int,
) -> int:
    """Sum chi(label,c) over c outside {0,excluded}."""

    if label == 0:
        return order - 2
    return -1 - character(label, excluded)


def finite_mgf_transform(order: int, label: int) -> float:
    """Walsh transform of the M_13 row amplitude f(r)."""

    return 0.0 if label == 0 else 1 / (order - 1)


def closed_record_three_spectrum(
    order: int,
    x: int,
    y: int,
    triple: tuple[int, int, int],
    mu: int,
    nu: int,
) -> float:
    """Closed twisted spectrum for horizontal pairs and a vertical triple.

    ``x`` and ``y`` are nonzero hidden-column differences.  ``triple`` is a
    canonical three-subset of the hidden-row group containing zero.  ``mu``
    and ``nu`` are full-cell Walsh frequencies.
    """

    if not (0 < x < order and 0 < y < order):
        raise ValueError((x, y))
    if triple[0] != 0 or len(set(triple)) != 3:
        raise ValueError(triple)
    mu_row, mu_column = divmod(mu, order)
    nu_row, nu_column = divmod(nu, order)
    c0 = character_sum_excluding_pair(order, mu_column, x)
    c1 = character_sum_excluding_pair(order, mu_column ^ y, x)

    # Translations whose vertical triple avoids row zero.  The two
    # chi(r,h) factors (one from M_35 and one from the Clifford twist) cancel.
    row_transform = sum(
        finite_mgf_transform(order, mu_row ^ value)
        for value in triple
    )
    if nu_row == 0:
        outside_sum = order - 3
    else:
        outside_sum = -sum(character(nu_row, value) for value in triple)
    # When the triple avoids row zero, the selected horizontal pair forms
    # one nontrivial even row of the quintic.  The record-three injection
    # formula contributes its conditional character average.
    even_row_c_sum = -(
        (1 + character(x, y)) * c0 + c1
    ) / (order - 3)
    nonincident = 0.0
    if nu_column == x:
        nonincident = (
            2 * order * even_row_c_sum * row_transform * outside_sum
        )

    # The three translations for which the triple contains row zero.  Its
    # column translation must then avoid the two selected-pair cells.
    column_label = nu_column ^ x
    if column_label == 0:
        z_sum = order - 2
    else:
        z_sum = -1 - character(column_label, y)
    incident = 0.0
    pair_phase = 1 + character(x, y)
    for h in triple:
        first = sum(
            finite_mgf_transform(order, mu_row ^ value)
            for value in triple
            if value != h
        )
        second = finite_mgf_transform(order, mu_row ^ h)
        incident += character(nu_row, h) * (
            pair_phase * c0 * first + 2 * c1 * second
        )
    denominator = order * (order - 1) * (order - 2)
    return (nonincident + z_sum * incident) / denominator


def frequency_class_representative(name: str) -> tuple[int, int]:
    representatives = {
        "00": (0, 0),
        "0*": (0, 1),
        "*0": (1, 0),
        "**0": (1, 2),
        "**1": (1, 1),
    }
    return representatives[name]


def frequency_class_multiplicity(order: int, name: str) -> int:
    if name == "00":
        return 1
    if name in ("0*", "*0"):
        return order - 1
    if name == "**0":
        return (order - 1) * (order // 2 - 1)
    if name == "**1":
        return (order - 1) * order // 2
    raise ValueError(name)


def closed_core_matrices(order: int) -> dict[tuple[str, str], np.ndarray]:
    """Build one exact core matrix for each of the 25 frequency orbits."""

    if order < 4 or order & (order - 1):
        raise ValueError(order)
    differences = tuple(range(1, order))
    triples = triple_orbit_representatives(order)
    names = ("00", "0*", "*0", "**0", "**1")
    result = {}
    for row_name in names:
        p_row, nu_row = frequency_class_representative(row_name)
        for column_name in names:
            p_column, nu_column = frequency_class_representative(column_name)
            nu = nu_row * order + nu_column
            matrix = np.empty(
                (len(differences) ** 2, len(triples)),
                dtype=float,
            )
            for triple_index, triple in enumerate(triples):
                triple_xor = xor_values(list(triple))
                mu = (p_row ^ triple_xor) * order + p_column
                for x_index, x in enumerate(differences):
                    for y_index, y in enumerate(differences):
                        matrix[x_index * len(differences) + y_index, triple_index] = (
                            closed_record_three_spectrum(
                                order,
                                x,
                                y,
                                triple,
                                mu,
                                nu,
                            )
                        )
            result[(row_name, column_name)] = matrix
    return result


def closed_symmetric_objective(
    order: int,
    cores: dict[tuple[str, str], np.ndarray],
    orthogonal_mass: float,
) -> float:
    differences = tuple(range(1, order))
    triples = triple_orbit_representatives(order)
    orthogonal_count = (order - 1) * (order // 2 - 1)
    nonorthogonal_count = (order - 1) * order // 2
    row_weights = np.asarray(
        [
            (
                orthogonal_mass / orthogonal_count
                if bilinear(x, y) == 0
                else (1 - orthogonal_mass) / nonorthogonal_count
            )
            for x in differences
            for y in differences
        ]
    )
    column_root = 1 / np.sqrt(len(triples))
    total = 0.0
    for key, core in cores.items():
        matrix = np.sqrt(row_weights)[:, None] * core * column_root
        nuclear = float(np.linalg.svd(matrix, compute_uv=False).sum())
        total += (
            frequency_class_multiplicity(order, key[0])
            * frequency_class_multiplicity(order, key[1])
            * nuclear
        )
    return total / order**4


def closed_symmetric_vertical_triple_witness(
    order: int,
) -> ClosedSymmetricVerticalTripleWitness:
    cores = closed_core_matrices(order)
    optimum = minimize_scalar(
        lambda mass: -closed_symmetric_objective(
            order,
            cores,
            float(mass),
        ),
        bounds=(1e-9, 1 - 1e-9),
        method="bounded",
        options={"xatol": 1e-12},
    )
    return ClosedSymmetricVerticalTripleWitness(
        order=order,
        orthogonal_mass=float(optimum.x),
        orthogonal_pairs=(order - 1) * (order // 2 - 1),
        nonorthogonal_pairs=(order - 1) * order // 2,
        triple_orbits=(order - 1) * (order - 2) // 6,
        frequency_classes=len(cores),
        coefficient=-float(optimum.fun),
    )


def build_data(order: int) -> SymmetricVerticalTripleData:
    if order not in (4, 8):
        raise ValueError("complete symmetric build is limited to q=4 or q=8")
    dimension = order**2
    differences = np.arange(1, order, dtype=np.int16)
    row_triples = triple_orbit_representatives(order)
    triples = tuple(
        tuple(value * order for value in triple)
        for triple in row_triples
    )
    triple_xors = np.asarray(
        [xor_values(list(triple)) for triple in triples],
        dtype=np.int16,
    )
    shape = (
        len(differences),
        len(differences),
        len(triples),
        dimension,
        dimension,
    )
    record_one_spectra = np.zeros(shape, dtype=float)
    record_three_spectra = np.zeros(shape, dtype=float)
    for x_index, x_value in enumerate(differences):
        x = int(x_value)
        middle = np.zeros(dimension, dtype=float)
        cubics: list[tuple[int, int, int] | None] = []
        records = np.zeros(dimension, dtype=np.int8)
        for s in range(dimension):
            if s in (0, x):
                cubics.append(None)
                continue
            cubic = tuple(sorted((0, x, s)))
            cubics.append(cubic)
            middle[s] = record_one_link_moment(order, (0,), cubic)
            records[s] = parity_record_size(order, cubic, axis=1)
        for y_index, y_value in enumerate(differences):
            y = int(y_value)
            for triple_index, triple in enumerate(triples):
                record_one = np.zeros((dimension, dimension), dtype=float)
                record_three = np.zeros_like(record_one)
                for s, cubic in enumerate(cubics):
                    if cubic is None or middle[s] == 0:
                        continue
                    for t in range(dimension):
                        shifted = tuple(value ^ t for value in triple)
                        if 0 in shifted or y in shifted:
                            continue
                        quintic = tuple(sorted((0, y) + shifted))
                        value = middle[s] * combined_link_moment(
                            order,
                            cubic,
                            quintic,
                        )
                        if records[s] == 1:
                            record_one[s, t] = value
                        elif records[s] == 3:
                            record_three[s, t] = value
                record_one_spectra[x_index, y_index, triple_index] = (
                    twisted_walsh(record_one)
                )
                record_three_spectra[x_index, y_index, triple_index] = (
                    twisted_walsh(record_three)
                )
    return SymmetricVerticalTripleData(
        order=order,
        differences=differences,
        triples=triples,
        triple_xors=triple_xors,
        record_one_spectra=record_one_spectra,
        record_three_spectra=record_three_spectra,
    )


def row_law(data: SymmetricVerticalTripleData, orthogonal_mass: float) -> np.ndarray:
    if not 0 <= orthogonal_mass <= 1:
        raise ValueError(orthogonal_mass)
    orthogonal = np.asarray(
        [
            [bilinear(int(x), int(y)) == 0 for y in data.differences]
            for x in data.differences
        ]
    )
    result = np.empty(orthogonal.shape, dtype=float)
    result[orthogonal] = orthogonal_mass / int(orthogonal.sum())
    result[~orthogonal] = (1 - orthogonal_mass) / int((~orthogonal).sum())
    return result


def evaluate_spectra(
    data: SymmetricVerticalTripleData,
    spectra: np.ndarray,
    orthogonal_mass: float,
) -> float:
    dimension = data.order**2
    row = row_law(data, orthogonal_mass)
    row_root = np.sqrt(row).reshape(-1)
    column_root = 1 / np.sqrt(len(data.triples))
    x_indices, y_indices = np.indices(row.shape)
    x_indices = x_indices.reshape(-1)
    y_indices = y_indices.reshape(-1)
    objective = 0.0
    for p in range(dimension):
        mu = np.bitwise_xor(p, data.triple_xors)
        for nu in range(dimension):
            core = spectra[
                x_indices[:, None],
                y_indices[:, None],
                np.arange(len(data.triples))[None, :],
                mu[None, :],
                nu,
            ]
            matrix = row_root[:, None] * core * column_root
            objective += float(np.linalg.svd(matrix, compute_uv=False).sum())
    return objective / dimension**2


def frequency_class_table(
    data: SymmetricVerticalTripleData,
    spectra: np.ndarray,
    orthogonal_mass: float,
) -> dict[tuple[str, str], tuple[int, float, float]]:
    """Group block nuclear norms by row- and column-frequency orbits."""

    order = data.order
    dimension = order**2
    row = row_law(data, orthogonal_mass)
    row_root = np.sqrt(row).reshape(-1)
    column_root = 1 / np.sqrt(len(data.triples))
    x_indices, y_indices = np.indices(row.shape)
    x_indices = x_indices.reshape(-1)
    y_indices = y_indices.reshape(-1)
    grouped: dict[tuple[str, str], list[float]] = {}
    for p in range(dimension):
        p_row, p_column = divmod(p, order)
        mu = np.bitwise_xor(p, data.triple_xors)
        for nu in range(dimension):
            nu_row, nu_column = divmod(nu, order)
            core = spectra[
                x_indices[:, None],
                y_indices[:, None],
                np.arange(len(data.triples))[None, :],
                mu[None, :],
                nu,
            ]
            value = float(
                np.linalg.svd(
                    row_root[:, None] * core * column_root,
                    compute_uv=False,
                ).sum()
            )
            key = (
                pair_frequency_class(p_row, nu_row),
                pair_frequency_class(p_column, nu_column),
            )
            grouped.setdefault(key, []).append(value)
    return {
        key: (len(values), min(values), max(values))
        for key, values in grouped.items()
    }


def symmetric_vertical_triple_witness(
    order: int,
) -> SymmetricVerticalTripleWitness:
    data = build_data(order)
    combined = data.record_one_spectra + data.record_three_spectra
    optimum = minimize_scalar(
        lambda mass: -evaluate_spectra(data, combined, float(mass)),
        bounds=(1e-8, 1 - 1e-8),
        method="bounded",
        options={"xatol": 1e-11},
    )
    mass = float(optimum.x)
    coefficient = -float(optimum.fun)
    record_one = evaluate_spectra(data, data.record_one_spectra, mass)
    record_three = evaluate_spectra(data, data.record_three_spectra, mass)
    row = row_law(data, mass)
    orthogonal = sum(
        bilinear(int(x), int(y)) == 0
        for x in data.differences
        for y in data.differences
    )
    return SymmetricVerticalTripleWitness(
        order=order,
        orthogonal_mass=mass,
        orthogonal_pairs=orthogonal,
        nonorthogonal_pairs=row.size - orthogonal,
        triple_orbits=len(data.triples),
        coefficient=coefficient,
        record_one_triangle=record_one,
        record_three_triangle=record_three,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--classes", action="store_true")
    parser.add_argument(
        "--closed",
        action="store_true",
        help="use the 25-class closed formula, including at q=16 or q=32",
    )
    parser.add_argument("orders", nargs="*", type=int, default=(4, 8))
    arguments = parser.parse_args()
    for order in arguments.orders:
        if arguments.closed:
            closed = closed_symmetric_vertical_triple_witness(order)
            print(
                "closed symmetric vertical-triple witness: "
                f"q={closed.order},N={closed.order**2},"
                f"orthogonal_mass={closed.orthogonal_mass:.12g},"
                f"pair_classes={closed.orthogonal_pairs}+"
                f"{closed.nonorthogonal_pairs},"
                f"triple_orbits={closed.triple_orbits},"
                f"frequency_classes={closed.frequency_classes},"
                f"coefficient={closed.coefficient:.12g}"
            )
            continue
        result = symmetric_vertical_triple_witness(order)
        print(
            "symmetric vertical-triple adjacent witness: "
            f"q={result.order},N={result.order**2},"
            f"orthogonal_mass={result.orthogonal_mass:.12g},"
            f"pair_classes={result.orthogonal_pairs}+"
            f"{result.nonorthogonal_pairs},"
            f"triple_orbits={result.triple_orbits},"
            f"record_one_triangle={result.record_one_triangle:.12g},"
            f"record_three_triangle={result.record_three_triangle:.12g},"
            f"coefficient={result.coefficient:.12g}"
        )
        if arguments.classes:
            data = build_data(order)
            table = frequency_class_table(
                data,
                data.record_one_spectra + data.record_three_spectra,
                result.orthogonal_mass,
            )
            for key, (count, minimum, maximum) in sorted(table.items()):
                print(
                    "frequency class: "
                    f"row={key[0]},column={key[1]},count={count},"
                    f"minimum={minimum:.12g},maximum={maximum:.12g},"
                    f"spread={maximum - minimum:.3g}"
                )


if __name__ == "__main__":
    main()
