#!/usr/bin/env python3
"""Complementary-row contraction for the reversed q64 middle pair.

For profile (1,1,5,3), the splits (0,1,3,1) and (0,1,2,2) share the same
complete-row bound.  Use the complete column as a transposed row for the
first and the declared row for the second.  Hadamard flatness cancels the
summed singleton, the exact fixed-pair quintic endpoint slice controls the
quintic completions, and only N-2 cubic singleton completions remain after
taking the universal middle-link maximum.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from math import sqrt
from pathlib import Path

from q64_paper_target_gate import (
    DIMENSION,
    MODES,
    ORDER,
    OptimizedLedger,
    ProfileSplit,
    optimize,
)
import occupation_compatible_sector_optimization as occupation
from middle_cubic_quintic_pair_contraction import middle_link_maxima
from q64_shifted_middle_pair_contraction import (
    inserted_coefficients as shifted_inserted_coefficients,
)


ROOT = Path(__file__).resolve().parents[1]
PROFILE = (1, 1, 5, 3)
SPLITS = ((0, 1, 3, 1), (0, 1, 2, 2))


@dataclass(frozen=True)
class ReversedMiddlePairContraction:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    quintic_fixed_pair_energy: float
    cubic_completions: int
    record_one_middle_maximum: float
    record_three_middle_maximum: float
    universal_middle_maximum: float
    row_energy_bound: float
    coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    reversed_pair_inserted: OptimizedLedger
    margin_improvement: float


def reversed_middle_pair_entries() -> tuple[ProfileSplit, ...]:
    reverse = tuple(reversed(PROFILE))
    entries: set[ProfileSplit] = set()
    for split in SPLITS:
        complement = tuple(
            degree - selected
            for degree, selected in zip(PROFILE, split, strict=True)
        )
        entries.update(
            {
                (PROFILE, split),
                (PROFILE, complement),
                (reverse, tuple(reversed(split))),
                (reverse, tuple(reversed(complement))),
            }
        )
    return tuple(sorted(entries))


def reversed_middle_pair_row_energy(order: int = ORDER) -> float:
    dimension = order * order
    quintic = occupation.endpoint_quintic_singleton_slice_energies(order)[2]
    middle = middle_link_maxima(order)[2]
    return (dimension - 2) * quintic * middle * middle


def reversed_middle_pair_coefficient(order: int = ORDER) -> float:
    return sqrt(reversed_middle_pair_row_energy(order))


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = shifted_inserted_coefficients()
    coefficient = reversed_middle_pair_coefficient()
    for entry in reversed_middle_pair_entries():
        result[entry] = coefficient
    return result


def diagnostic() -> ReversedMiddlePairContraction:
    entries = reversed_middle_pair_entries()
    quintic = occupation.endpoint_quintic_singleton_slice_energies(ORDER)[2]
    record_one, record_three, middle = middle_link_maxima(ORDER)
    row_energy = reversed_middle_pair_row_energy()
    previous = optimize(mapped_coefficients=shifted_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    previous_proved = 224
    return ReversedMiddlePairContraction(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        quintic_fixed_pair_energy=quintic,
        cubic_completions=DIMENSION - 2,
        record_one_middle_maximum=record_one,
        record_three_middle_maximum=record_three,
        universal_middle_maximum=middle,
        row_energy_bound=row_energy,
        coefficient=sqrt(row_energy),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        reversed_pair_inserted=inserted,
        margin_improvement=(
            inserted.margin_to_one_third - previous.margin_to_one_third
        ),
    )


def artifact_text(result: ReversedMiddlePairContraction) -> str:
    payload = {
        "schema": "round4_q64_reversed_middle_pair_contraction_v1",
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal complementary complete-row Schur-feature "
            "theorem using Hadamard flatness, the exact quintic fixed-pair "
            "slice, and the universal middle-link maximum; floating q64 "
            "Perron insertion; one batch only"
        ),
        "remaining_open": (
            "656 balanced entries plus interval certification and the adaptive lift"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic reversed middle-pair contraction",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_reversed_middle_pair_contraction.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 reversed middle-pair contraction: "
        f"entries={result.closed_entries},"
        f"coefficient={result.coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.reversed_pair_inserted.total:.12g},"
        f"margin={result.reversed_pair_inserted.margin_to_one_third:.12g},"
        f"margin_gain={result.margin_improvement:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
