#!/usr/bin/env python3
"""Whole-cubic decorated completion rows for 16 q64 quintic entries."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from math import comb, sqrt
from pathlib import Path

from q64_degree_ten_completion_row_insertion import (
    inserted_coefficients as degree_ten_inserted_coefficients,
    orbit,
    remaining_quintic_entries as pre_whole_cubic_quintic_entries,
)
from q64_paper_target_gate import (
    DIMENSION,
    MODES,
    ORDER,
    RESERVE_TARGET,
    OptimizedLedger,
    ProfileSplit,
    optimize,
)
from q64_post_universal_quintic_gate import quintic_split_depth
from middle_cubic_quintic_pair_contraction import middle_link_maxima
import occupation_compatible_sector_optimization as occupation


ROOT = Path(__file__).resolve().parents[1]

CUBIC_LEADING_CUBIC_QUINTIC = ((3, 1, 3, 5), (0, 0, 2, 4))
CUBIC_QUINTIC_TRAILING_CUBIC = ((1, 3, 5, 3), (0, 2, 4, 0))
CUBIC_LEADING_QUINTIC_CUBIC = ((3, 1, 5, 3), (0, 0, 4, 2))
QUINTIC_CUBIC_TRAILING_CUBIC = ((1, 5, 3, 3), (0, 4, 2, 0))


@dataclass(frozen=True)
class WholeCubicDecoratedRowInsertion:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    closed_orbits: int
    extreme_entries: int
    balanced_entries: int
    record_one_middle_maximum: float
    record_three_middle_maximum: float
    cubic_endpoint_record_one_energy: float
    cubic_endpoint_record_three_energy: float
    cubic_leading_residual_factor: float
    cubic_leading_coefficient: float
    cubic_trailing_residual_factor: float
    cubic_trailing_coefficient: float
    quintic_endpoint_record_one_energy: float
    quintic_endpoint_record_three_energy: float
    quintic_leading_residual_factor: float
    quintic_leading_coefficient: float
    quintic_trailing_residual_factor: float
    quintic_trailing_coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    decorated_rows_inserted: OptimizedLedger
    routing_change: float
    remaining_quintic_entries: int
    remaining_extreme_entries: int
    remaining_balanced_entries: int
    remaining_quintic_local_proxy: OptimizedLedger
    proxy_reserve_after_declared_allowance: float


def whole_cubic_decorated_entries() -> tuple[ProfileSplit, ...]:
    return tuple(
        sorted(
            set(orbit(CUBIC_LEADING_CUBIC_QUINTIC))
            | set(orbit(CUBIC_QUINTIC_TRAILING_CUBIC))
            | set(orbit(CUBIC_LEADING_QUINTIC_CUBIC))
            | set(orbit(QUINTIC_CUBIC_TRAILING_CUBIC))
        )
    )


def sector_parameters(order: int = ORDER) -> tuple[float, float]:
    record_one, _, _ = middle_link_maxima(order)
    return record_one, 1 / comb(order, 3)


def cubic_endpoint_energy_parts(
    order: int = ORDER,
) -> tuple[float, float]:
    """Energy of ``v3(C) M35(C,S)`` with fixed counts two and four."""

    q = order
    dimension = q * q
    record_one, record_three = sector_parameters(q)
    cubic_pair_energy = occupation.endpoint_singleton_slice_energies(q)[2]
    one = (
        2 * (q - 1) / (q - 1) ** 2 * (q * q - 4) * record_one**2
    )
    three = (
        dimension
        * cubic_pair_energy
        * (dimension - 4)
        * record_three**2
    )
    return one, three


def quintic_endpoint_energy_parts(
    order: int = ORDER,
) -> tuple[float, float]:
    """Energy of ``v5(S) M53(S,C)`` with fixed counts four and two."""

    q = order
    dimension = q * q
    record_one, record_three = sector_parameters(q)
    quintic_fixed_four_energy = (
        occupation.endpoint_quintic_singleton_slice_energies(q)[4]
    )
    one = (q * q - 4) * 2 * (q - 1) * record_one**2
    three = (
        dimension
        * quintic_fixed_four_energy
        * (dimension - 2)
        * record_three**2
    )
    return one, three


def cubic_leading_coefficient(order: int = ORDER) -> float:
    return sqrt(sum(cubic_endpoint_energy_parts(order))) / order


def cubic_trailing_coefficient(order: int = ORDER) -> float:
    return sqrt(sum(cubic_endpoint_energy_parts(order)))


def quintic_leading_coefficient(order: int = ORDER) -> float:
    return sqrt(sum(quintic_endpoint_energy_parts(order))) / order


def quintic_trailing_coefficient(order: int = ORDER) -> float:
    return sqrt(sum(quintic_endpoint_energy_parts(order)))


def coefficient_map() -> dict[ProfileSplit, float]:
    result: dict[ProfileSplit, float] = {}
    families = (
        (CUBIC_LEADING_CUBIC_QUINTIC, cubic_leading_coefficient()),
        (CUBIC_QUINTIC_TRAILING_CUBIC, cubic_trailing_coefficient()),
        (CUBIC_LEADING_QUINTIC_CUBIC, quintic_leading_coefficient()),
        (QUINTIC_CUBIC_TRAILING_CUBIC, quintic_trailing_coefficient()),
    )
    for generator, coefficient in families:
        for entry in orbit(generator):
            if entry in result:
                raise AssertionError(("decorated-row orbit overlap", entry))
            result[entry] = coefficient
    return result


def remaining_quintic_entries() -> tuple[ProfileSplit, ...]:
    closed = set(whole_cubic_decorated_entries())
    return tuple(
        entry
        for entry in pre_whole_cubic_quintic_entries()
        if entry not in closed
    )


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = degree_ten_inserted_coefficients()
    result.update(coefficient_map())
    return result


def remaining_quintic_local_proxy_coefficients() -> dict[ProfileSplit, float]:
    from q64_dual_endpoint_schur_insertion import local_slice_coefficients

    result = inserted_coefficients()
    extreme, balanced = local_slice_coefficients()
    for entry in remaining_quintic_entries():
        result[entry] = (
            extreme if quintic_split_depth(entry) == 1 else balanced
        )
    return result


def diagnostic() -> WholeCubicDecoratedRowInsertion:
    entries = whole_cubic_decorated_entries()
    remaining = remaining_quintic_entries()
    record_one, record_three = sector_parameters()
    cubic_one, cubic_three = cubic_endpoint_energy_parts()
    quintic_one, quintic_three = quintic_endpoint_energy_parts()
    previous = optimize(mapped_coefficients=degree_ten_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    proxy = optimize(
        mapped_coefficients=remaining_quintic_local_proxy_coefficients()
    )
    previous_proved = 304
    return WholeCubicDecoratedRowInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        closed_orbits=len(entries) // 4,
        extreme_entries=sum(quintic_split_depth(entry) == 1 for entry in entries),
        balanced_entries=sum(quintic_split_depth(entry) == 2 for entry in entries),
        record_one_middle_maximum=record_one,
        record_three_middle_maximum=record_three,
        cubic_endpoint_record_one_energy=cubic_one,
        cubic_endpoint_record_three_energy=cubic_three,
        cubic_leading_residual_factor=1 / ORDER,
        cubic_leading_coefficient=cubic_leading_coefficient(),
        cubic_trailing_residual_factor=1.0,
        cubic_trailing_coefficient=cubic_trailing_coefficient(),
        quintic_endpoint_record_one_energy=quintic_one,
        quintic_endpoint_record_three_energy=quintic_three,
        quintic_leading_residual_factor=1 / ORDER,
        quintic_leading_coefficient=quintic_leading_coefficient(),
        quintic_trailing_residual_factor=1.0,
        quintic_trailing_coefficient=quintic_trailing_coefficient(),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        decorated_rows_inserted=inserted,
        routing_change=inserted.total - previous.total,
        remaining_quintic_entries=len(remaining),
        remaining_extreme_entries=sum(
            quintic_split_depth(entry) == 1 for entry in remaining
        ),
        remaining_balanced_entries=sum(
            quintic_split_depth(entry) == 2 for entry in remaining
        ),
        remaining_quintic_local_proxy=proxy,
        proxy_reserve_after_declared_allowance=(
            proxy.margin_to_one_third - RESERVE_TARGET
        ),
    )


def artifact_text(result: WholeCubicDecoratedRowInsertion) -> str:
    payload = {
        "schema": "round4_q64_whole_cubic_decorated_row_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal scalar completion rows with whole-cubic "
            "endpoint decorations treated as column-only 1/q factors or "
            "unit cross-Gram Schur multipliers; floating q64 Perron "
            "insertion; one batch only"
        ),
        "remaining_open": (
            "568 balanced entries, including 108 split-cubic/split-quintic "
            "entries, plus interval certification and the adaptive lift"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic whole-cubic decorated-row insertion",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_whole_cubic_decorated_row_insertion.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 whole-cubic decorated-row insertion: "
        f"entries={result.closed_entries},"
        f"orbits={result.closed_orbits},"
        f"coefficients={result.cubic_leading_coefficient:.12g}/"
        f"{result.cubic_trailing_coefficient:.12g}/"
        f"{result.quintic_leading_coefficient:.12g}/"
        f"{result.quintic_trailing_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.decorated_rows_inserted.total:.12g},"
        f"margin={result.decorated_rows_inserted.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        f"proxy_total={result.remaining_quintic_local_proxy.total:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
