#!/usr/bin/env python3
"""Fixed-singleton complementary-row contraction at q=64.

For profile (1,1,5,3) and split (0,0,3,2), use the complementary split.
The row then fixes both singleton blocks, a quintic pair, and one cubic
cell.  Hadamard flatness contributes 1/N, the exact quintic fixed-pair
slice controls its triple completions, and C(N-1,2) cubic pair completions
remain after the universal middle-link maximum is taken.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from math import comb, sqrt
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
from q64_reversed_middle_pair_contraction import (
    inserted_coefficients as reversed_inserted_coefficients,
)


ROOT = Path(__file__).resolve().parents[1]
PROFILE = (1, 1, 5, 3)
SPLIT = (0, 0, 3, 2)


@dataclass(frozen=True)
class FixedSingletonPairContraction:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    quintic_fixed_pair_energy: float
    cubic_pair_completions: int
    hadamard_squared_factor: float
    universal_middle_maximum: float
    row_energy_bound: float
    coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    fixed_singleton_inserted: OptimizedLedger
    margin_improvement: float


def fixed_singleton_pair_entries() -> tuple[ProfileSplit, ...]:
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


def fixed_singleton_pair_row_energy(order: int = ORDER) -> float:
    dimension = order * order
    quintic = occupation.endpoint_quintic_singleton_slice_energies(order)[2]
    middle = middle_link_maxima(order)[2]
    return quintic * comb(dimension - 1, 2) * middle * middle / dimension


def fixed_singleton_pair_coefficient(order: int = ORDER) -> float:
    return sqrt(fixed_singleton_pair_row_energy(order))


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = reversed_inserted_coefficients()
    coefficient = fixed_singleton_pair_coefficient()
    for entry in fixed_singleton_pair_entries():
        result[entry] = coefficient
    return result


def diagnostic() -> FixedSingletonPairContraction:
    entries = fixed_singleton_pair_entries()
    quintic = occupation.endpoint_quintic_singleton_slice_energies(ORDER)[2]
    middle = middle_link_maxima(ORDER)[2]
    row_energy = fixed_singleton_pair_row_energy()
    previous = optimize(mapped_coefficients=reversed_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    previous_proved = 232
    return FixedSingletonPairContraction(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        quintic_fixed_pair_energy=quintic,
        cubic_pair_completions=comb(DIMENSION - 1, 2),
        hadamard_squared_factor=1 / DIMENSION,
        universal_middle_maximum=middle,
        row_energy_bound=row_energy,
        coefficient=sqrt(row_energy),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        fixed_singleton_inserted=inserted,
        margin_improvement=(
            inserted.margin_to_one_third - previous.margin_to_one_third
        ),
    )


def artifact_text(result: FixedSingletonPairContraction) -> str:
    payload = {
        "schema": "round4_q64_fixed_singleton_pair_contraction_v1",
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal complementary complete-row Schur-feature "
            "theorem using Hadamard flatness, the exact quintic fixed-pair "
            "slice, a counted cubic completion, and the universal middle-link "
            "maximum; floating q64 Perron insertion; one batch only"
        ),
        "remaining_open": (
            "652 balanced entries plus interval certification and the adaptive lift"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic fixed-singleton pair contraction",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_fixed_singleton_pair_contraction.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 fixed-singleton pair contraction: "
        f"entries={result.closed_entries},"
        f"coefficient={result.coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.fixed_singleton_inserted.total:.12g},"
        f"margin={result.fixed_singleton_inserted.margin_to_one_third:.12g},"
        f"margin_gain={result.margin_improvement:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
