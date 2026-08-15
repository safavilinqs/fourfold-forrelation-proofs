#!/usr/bin/env python3
"""Valid non-invariant two-axis screen for the leading physical orbit.

The invariant opposite-endpoint witness uses two vertical pair differences
and a vertical triple orbit. The unrestricted ``q=4`` optimizer suggests
testing a small horizontal component in the quintic selected-pair
difference. This module evaluates exactly that physical family at ``q=32``.

The law remains uniform under independent row-label symmetries inside each
chosen axis, so the full frequency sum reduces to 25 representative blocks.
Unlike the fourteen-class invariant reduction, this reduction is also
checked against the unreduced ``q=4`` formula for a genuinely non-invariant
law.

The result is a numerical physical lower witness, not an arbitrary-law upper
bound or an interval certificate.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from functools import lru_cache
from json import dumps
from math import sqrt
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent / "open_problem_forr4_passive_floor_consolidation_round_3" / "searches"
)
sys.path.insert(0, str(ROUND3_SEARCHES))

from opposite_endpoint_orbit_scan import (  # noqa: E402
    cubic_response,
    quintic_response,
    support_xor,
    triple_orbit_representatives,
    walsh_transform,
)


ORDER = 32
PHYSICAL_GATE = 0.0414623182965146


@dataclass(frozen=True)
class TwoAxisData:
    order: int
    differences: np.ndarray
    vertical_triples: tuple[tuple[int, int, int], ...]
    cubic_spectrum: np.ndarray
    quintic_spectrum: np.ndarray


@dataclass(frozen=True)
class TwoAxisSearch:
    order: int
    invariant_equal_mass: float
    invariant_coefficient: float
    best_horizontal_mass: float
    best_coefficient: float
    absolute_improvement: float
    relative_improvement: float
    physical_gate: float
    gate_headroom: float
    decision: str
    evidence_label: str


def character_value(left: int, right: int) -> int:
    return -1 if (left & right).bit_count() % 2 else 1


@lru_cache(maxsize=None)
def build_data(order: int) -> TwoAxisData:
    if order < 4 or order & (order - 1):
        raise ValueError(("power-of-two order at least four required", order))
    q = order
    dimension = q * q
    vertical = tuple(row * q for row in range(1, q))
    horizontal = tuple(range(1, q))
    differences = np.asarray(vertical + horizontal, dtype=np.int64)
    vertical_triples = tuple(
        tuple(row * q for row in triple) for triple in triple_orbit_representatives(q)
    )
    cubic = np.asarray(
        [
            walsh_transform(cubic_response(q, int(difference), False))
            for difference in differences
        ]
    )
    quintic = np.empty(
        (len(differences), len(vertical_triples), dimension), dtype=float
    )
    frequencies = tuple(range(dimension))
    for difference_index, difference in enumerate(differences):
        for triple_index, triple in enumerate(vertical_triples):
            triple_xor = support_xor(triple)
            characters = np.asarray(
                [character_value(triple_xor, frequency) for frequency in frequencies]
            )
            quintic[difference_index, triple_index] = (
                walsh_transform(quintic_response(q, int(difference), triple, False))
                * characters
            )
    return TwoAxisData(
        order=q,
        differences=differences,
        vertical_triples=vertical_triples,
        cubic_spectrum=cubic,
        quintic_spectrum=quintic,
    )


def frequency_classes(order: int) -> tuple[tuple[int, int, int], ...]:
    q = order
    return (
        (0, 0, 1),
        (0, 1, q - 1),
        (1, 0, q - 1),
        (1, 2, (q - 1) * (q // 2 - 1)),
        (1, 1, (q - 1) * (q // 2)),
    )


def row_law(
    order: int,
    equal_vertical_mass: float,
    horizontal_quintic_mass: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not 0 <= equal_vertical_mass <= 1:
        raise ValueError(("equal vertical mass", equal_vertical_mass))
    if not 0 <= horizontal_quintic_mass <= 1:
        raise ValueError(("horizontal mass", horizontal_quintic_mass))
    count = order - 1
    vertical_pair_count = count * count
    vertical_indices = np.repeat(np.arange(count), count)
    vertical_quintic = np.tile(np.arange(count), count)
    equal = vertical_indices == vertical_quintic
    vertical_probabilities = np.where(
        equal,
        equal_vertical_mass / count,
        (1 - equal_vertical_mass) / (count * (count - 1)),
    ) * (1 - horizontal_quintic_mass)

    horizontal_indices = np.repeat(np.arange(count), count)
    horizontal_quintic = np.tile(np.arange(count, 2 * count), count)
    horizontal_probabilities = np.full(
        vertical_pair_count,
        horizontal_quintic_mass / vertical_pair_count,
    )
    cubic_indices = np.concatenate((vertical_indices, horizontal_indices))
    quintic_indices = np.concatenate((vertical_quintic, horizontal_quintic))
    probabilities = np.concatenate((vertical_probabilities, horizontal_probabilities))
    active = probabilities > 0
    return (
        cubic_indices[active],
        quintic_indices[active],
        probabilities[active],
    )


def coefficient(
    order: int,
    equal_vertical_mass: float,
    horizontal_quintic_mass: float,
) -> float:
    data = build_data(order)
    q = order
    dimension = q * q
    cubic_indices, quintic_indices, probabilities = row_law(
        q, equal_vertical_mass, horizontal_quintic_mass
    )
    differences = data.differences[cubic_indices]
    row_amplitudes = np.sqrt(probabilities)
    column_amplitude = 1 / sqrt(len(data.vertical_triples))
    columns = np.arange(len(data.vertical_triples))
    total = 0.0
    for alpha_row, gamma_row, row_multiplicity in frequency_classes(q):
        for alpha_column, gamma_column, column_multiplicity in frequency_classes(q):
            alpha = alpha_row * q + alpha_column
            gamma = gamma_row * q + gamma_column
            shifted = np.bitwise_xor(gamma, differences)
            endpoint = data.quintic_spectrum[
                quintic_indices[:, None],
                columns[None, :],
                shifted[:, None],
            ]
            matrix = (
                (row_amplitudes * data.cubic_spectrum[cubic_indices, alpha])[:, None]
                * endpoint
                * column_amplitude
            )
            total += (
                row_multiplicity
                * column_multiplicity
                * float(np.linalg.svd(matrix, compute_uv=False).sum())
            )
    return total / dimension**3


def search(order: int = ORDER) -> TwoAxisSearch:
    invariant_mass_search = minimize_scalar(
        lambda mass: -coefficient(order, float(mass), 0.0),
        bounds=(0.0, 0.1),
        method="bounded",
        options={"xatol": 1e-12},
    )
    invariant_mass = float(invariant_mass_search.x)
    invariant = -float(invariant_mass_search.fun)
    horizontal_search = minimize_scalar(
        lambda mass: -coefficient(order, invariant_mass, float(mass)),
        bounds=(0.0, 0.01),
        method="bounded",
        options={"xatol": 1e-12},
    )
    horizontal_mass = float(horizontal_search.x)
    best = -float(horizontal_search.fun)
    return TwoAxisSearch(
        order=order,
        invariant_equal_mass=invariant_mass,
        invariant_coefficient=invariant,
        best_horizontal_mass=horizontal_mass,
        best_coefficient=best,
        absolute_improvement=best - invariant,
        relative_improvement=best / invariant - 1,
        physical_gate=PHYSICAL_GATE,
        gate_headroom=PHYSICAL_GATE - best,
        decision="does_not_cross_leading_orbit_kill_gate",
        evidence_label=(
            "valid physical lower-witness family with exact response formulas "
            "and floating nuclear norms; not interval certified"
        ),
    )


def artifact_text(result: TwoAxisSearch) -> str:
    return (
        dumps(
            {
                "schema": "round4_opposite_endpoint_two_axis_screen_v1",
                **asdict(result),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic two-axis screen artifact",
    )
    arguments = parser.parse_args()
    result = search()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "opposite_endpoint_two_axis_screen.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "opposite-endpoint two-axis screen: "
        f"q={result.order},"
        f"invariant_mass={result.invariant_equal_mass:.15g},"
        f"invariant={result.invariant_coefficient:.15g},"
        f"horizontal_mass={result.best_horizontal_mass:.15g},"
        f"best={result.best_coefficient:.15g},"
        f"improvement={result.absolute_improvement:.15g},"
        f"relative_improvement={result.relative_improvement:.15g},"
        f"gate={result.physical_gate:.15g},"
        f"headroom={result.gate_headroom:.15g},"
        f"decision={result.decision}"
    )


if __name__ == "__main__":
    main()
