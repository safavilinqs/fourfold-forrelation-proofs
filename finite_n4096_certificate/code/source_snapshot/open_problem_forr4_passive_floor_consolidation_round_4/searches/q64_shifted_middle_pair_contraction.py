#!/usr/bin/env python3
"""Shifted middle cubic-quintic row contraction at q=64.

This is the 1|2 cubic / 3|2 quintic shift of the accepted Round 3 middle
pair theorem.  A complete row fixes one cubic cell, three quintic cells,
and the final singleton.  Taking the universal middle-link maximum separates
the exact endpoint slice sums, and the complete row becomes a Schur feature.
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
from q64_universal_multicubic_insertion import (
    inserted_coefficients as universal_inserted_coefficients,
)


ROOT = Path(__file__).resolve().parents[1]
PROFILE = (1, 3, 5, 1)
SPLIT = (0, 1, 3, 1)


@dataclass(frozen=True)
class ShiftedMiddlePairContraction:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    cubic_fixed_singleton_energy: float
    quintic_fixed_triple_energy: float
    record_one_middle_maximum: float
    record_three_middle_maximum: float
    universal_middle_maximum: float
    row_energy_bound: float
    coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    shifted_pair_inserted: OptimizedLedger
    margin_improvement: float


def shifted_middle_pair_entries() -> tuple[ProfileSplit, ...]:
    complement = tuple(
        degree - selected
        for degree, selected in zip(PROFILE, SPLIT, strict=True)
    )
    reverse = tuple(reversed(PROFILE))
    return tuple(
        sorted(
            {
                (PROFILE, SPLIT),
                (PROFILE, complement),
                (reverse, tuple(reversed(SPLIT))),
                (reverse, tuple(reversed(complement))),
            }
        )
    )


def shifted_middle_pair_row_energy(order: int = ORDER) -> float:
    dimension = order * order
    cubic = occupation.endpoint_singleton_slice_energies(order)[1]
    quintic = occupation.endpoint_quintic_singleton_slice_energies(order)[3]
    middle = middle_link_maxima(order)[2]
    return dimension * cubic * quintic * middle * middle


def shifted_middle_pair_coefficient(order: int = ORDER) -> float:
    return sqrt(shifted_middle_pair_row_energy(order))


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = universal_inserted_coefficients()
    coefficient = shifted_middle_pair_coefficient()
    for entry in shifted_middle_pair_entries():
        result[entry] = coefficient
    return result


def diagnostic() -> ShiftedMiddlePairContraction:
    entries = shifted_middle_pair_entries()
    cubic = occupation.endpoint_singleton_slice_energies(ORDER)[1]
    quintic = occupation.endpoint_quintic_singleton_slice_energies(ORDER)[3]
    record_one, record_three, middle = middle_link_maxima(ORDER)
    row_energy = shifted_middle_pair_row_energy()
    previous = optimize(mapped_coefficients=universal_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    previous_proved = 220
    return ShiftedMiddlePairContraction(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        cubic_fixed_singleton_energy=cubic,
        quintic_fixed_triple_energy=quintic,
        record_one_middle_maximum=record_one,
        record_three_middle_maximum=record_three,
        universal_middle_maximum=middle,
        row_energy_bound=row_energy,
        coefficient=sqrt(row_energy),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        shifted_pair_inserted=inserted,
        margin_improvement=(
            inserted.margin_to_one_third - previous.margin_to_one_third
        ),
    )


def artifact_text(result: ShiftedMiddlePairContraction) -> str:
    payload = {
        "schema": "round4_q64_shifted_middle_pair_contraction_v1",
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal complete-row Schur-feature theorem using "
            "exact cubic/quintic endpoint slices and the universal middle-link "
            "maximum; floating q64 Perron insertion; one batch only"
        ),
        "remaining_open": (
            "664 balanced entries plus interval certification and the adaptive lift"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic shifted middle-pair contraction",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_shifted_middle_pair_contraction.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 shifted middle-pair contraction: "
        f"entries={result.closed_entries},"
        f"coefficient={result.coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.shifted_pair_inserted.total:.12g},"
        f"margin={result.shifted_pair_inserted.margin_to_one_third:.12g},"
        f"margin_gain={result.margin_improvement:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
