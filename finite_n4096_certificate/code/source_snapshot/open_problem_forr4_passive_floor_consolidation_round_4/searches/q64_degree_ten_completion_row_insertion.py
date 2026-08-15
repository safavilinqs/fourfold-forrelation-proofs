#!/usr/bin/env python3
"""Completion-row contractions for 12 remaining q64 degree-ten entries."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from math import comb, sqrt
from pathlib import Path

from q64_decorated_adjacent_row_insertion import (
    inserted_coefficients as decorated_inserted_coefficients,
    remaining_quintic_entries as pre_degree_ten_quintic_entries,
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

LEFT_DOUBLE_SINGLETON = ((1, 1, 3, 5), (0, 0, 2, 3))
REVERSED_DOUBLE_SINGLETON = ((1, 1, 5, 3), (0, 0, 4, 1))
DOUBLE_ENDPOINT_ONE_FOUR = ((1, 3, 5, 1), (0, 1, 4, 0))


@dataclass(frozen=True)
class DegreeTenCompletionRowInsertion:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    closed_orbits: int
    extreme_entries: int
    balanced_entries: int
    record_one_middle_maximum: float
    record_three_middle_maximum: float
    left_double_singleton_record_one_energy: float
    left_double_singleton_record_three_energy: float
    left_double_singleton_residual_factor: float
    left_double_singleton_coefficient: float
    reversed_double_singleton_record_one_energy: float
    reversed_double_singleton_record_three_energy: float
    reversed_double_singleton_residual_factor: float
    reversed_double_singleton_coefficient: float
    double_endpoint_record_one_energy: float
    double_endpoint_record_three_energy: float
    double_endpoint_coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    degree_ten_inserted: OptimizedLedger
    routing_margin_improvement: float
    remaining_quintic_entries: int
    remaining_extreme_entries: int
    remaining_balanced_entries: int
    remaining_quintic_local_proxy: OptimizedLedger
    proxy_reserve_after_declared_allowance: float


def orbit(entry: ProfileSplit) -> tuple[ProfileSplit, ...]:
    profile, split = entry
    complement = tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )
    reverse = tuple(reversed(profile))
    return tuple(
        sorted(
            {
                (profile, split),
                (profile, complement),
                (reverse, tuple(reversed(split))),
                (reverse, tuple(reversed(complement))),
            }
        )
    )


def left_double_singleton_entries() -> tuple[ProfileSplit, ...]:
    return orbit(LEFT_DOUBLE_SINGLETON)


def reversed_double_singleton_entries() -> tuple[ProfileSplit, ...]:
    return orbit(REVERSED_DOUBLE_SINGLETON)


def double_endpoint_one_four_entries() -> tuple[ProfileSplit, ...]:
    return orbit(DOUBLE_ENDPOINT_ONE_FOUR)


def degree_ten_completion_entries() -> tuple[ProfileSplit, ...]:
    return tuple(
        sorted(
            set(left_double_singleton_entries())
            | set(reversed_double_singleton_entries())
            | set(double_endpoint_one_four_entries())
        )
    )


def sector_parameters(order: int = ORDER) -> tuple[float, float]:
    record_one, _, _ = middle_link_maxima(order)
    record_three = 1 / comb(order, 3)
    return record_one, record_three


def left_double_singleton_energy_parts(
    order: int = ORDER,
) -> tuple[float, float]:
    """Scalar-row energies for ``(1,1,3,5):(0,0,2,3)``.

    The record-one term is the accepted fixed-pair/triple incidence bound.
    In record three, the exact cubic fixed-pair endpoint slice is combined
    with the total number of quintic pair completions.  The residual product
    of the two singleton Walsh links contributes ``1/q`` after the row norm.
    """

    q = order
    dimension = q * q
    record_one, record_three = sector_parameters(q)
    cubic_pair_incidence = 2 * (q - 1)
    quintic_triple_incidence = occupation.middle_quintic_incidence_bound(q, 3)
    cubic_pair_energy = occupation.endpoint_singleton_slice_energies(q)[2]
    one = (
        cubic_pair_incidence
        / (q - 1) ** 2
        * quintic_triple_incidence
        * record_one**2
    )
    three = (
        dimension
        * cubic_pair_energy
        * comb(dimension - 3, 2)
        * record_three**2
    )
    return one, three


def left_double_singleton_coefficient(order: int = ORDER) -> float:
    return sqrt(sum(left_double_singleton_energy_parts(order))) / order


def reversed_double_singleton_energy_parts(
    order: int = ORDER,
) -> tuple[float, float]:
    """Scalar-row energies for ``(1,1,5,3):(0,0,4,1)``."""

    q = order
    dimension = q * q
    record_one, record_three = sector_parameters(q)
    quintic_fixed_four_incidence = q * q - 4
    cubic_fixed_one_incidence = 3 * (q - 1) ** 2
    one = (
        quintic_fixed_four_incidence
        * cubic_fixed_one_incidence
        * record_one**2
    )
    three = (
        (dimension - 4)
        * comb(dimension - 1, 2)
        * record_three**2
    )
    return one, three


def reversed_double_singleton_coefficient(order: int = ORDER) -> float:
    return sqrt(sum(reversed_double_singleton_energy_parts(order))) / order


def double_endpoint_one_four_energy_parts(
    order: int = ORDER,
) -> tuple[float, float]:
    """Scalar-row energies for ``(1,3,5,1):(0,1,4,0)``."""

    q = order
    dimension = q * q
    record_one, record_three = sector_parameters(q)
    cubic_fixed_one_energy = occupation.endpoint_singleton_slice_energies(q)[1]
    quintic_fixed_four_energy = (
        occupation.endpoint_quintic_singleton_slice_energies(q)[4]
    )
    one = 3 * (q * q - 4) * record_one**2
    three = (
        dimension**2
        * cubic_fixed_one_energy
        * quintic_fixed_four_energy
        * record_three**2
    )
    return one, three


def double_endpoint_one_four_coefficient(order: int = ORDER) -> float:
    return sqrt(sum(double_endpoint_one_four_energy_parts(order)))


def coefficient_map() -> dict[ProfileSplit, float]:
    result: dict[ProfileSplit, float] = {}
    for entry in left_double_singleton_entries():
        result[entry] = left_double_singleton_coefficient()
    for entry in reversed_double_singleton_entries():
        result[entry] = reversed_double_singleton_coefficient()
    for entry in double_endpoint_one_four_entries():
        result[entry] = double_endpoint_one_four_coefficient()
    return result


def remaining_quintic_entries() -> tuple[ProfileSplit, ...]:
    closed = set(degree_ten_completion_entries())
    return tuple(
        entry
        for entry in pre_degree_ten_quintic_entries()
        if entry not in closed
    )


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = decorated_inserted_coefficients()
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


def diagnostic() -> DegreeTenCompletionRowInsertion:
    entries = degree_ten_completion_entries()
    remaining = remaining_quintic_entries()
    record_one, record_three = sector_parameters()
    left_one, left_three = left_double_singleton_energy_parts()
    reversed_one, reversed_three = reversed_double_singleton_energy_parts()
    endpoint_one, endpoint_three = double_endpoint_one_four_energy_parts()
    previous = optimize(mapped_coefficients=decorated_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    proxy = optimize(
        mapped_coefficients=remaining_quintic_local_proxy_coefficients()
    )
    previous_proved = 292
    return DegreeTenCompletionRowInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        closed_orbits=len(entries) // 4,
        extreme_entries=sum(quintic_split_depth(entry) == 1 for entry in entries),
        balanced_entries=sum(quintic_split_depth(entry) == 2 for entry in entries),
        record_one_middle_maximum=record_one,
        record_three_middle_maximum=record_three,
        left_double_singleton_record_one_energy=left_one,
        left_double_singleton_record_three_energy=left_three,
        left_double_singleton_residual_factor=1 / ORDER,
        left_double_singleton_coefficient=left_double_singleton_coefficient(),
        reversed_double_singleton_record_one_energy=reversed_one,
        reversed_double_singleton_record_three_energy=reversed_three,
        reversed_double_singleton_residual_factor=1 / ORDER,
        reversed_double_singleton_coefficient=(
            reversed_double_singleton_coefficient()
        ),
        double_endpoint_record_one_energy=endpoint_one,
        double_endpoint_record_three_energy=endpoint_three,
        double_endpoint_coefficient=double_endpoint_one_four_coefficient(),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        degree_ten_inserted=inserted,
        routing_margin_improvement=(
            inserted.margin_to_one_third - previous.margin_to_one_third
        ),
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


def artifact_text(result: DegreeTenCompletionRowInsertion) -> str:
    payload = {
        "schema": "round4_q64_degree_ten_completion_row_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal scalar completion-row contractions using "
            "record-one incidence, record-three completion counts, exact "
            "endpoint slices, and unitary Walsh residuals; floating q64 "
            "Perron insertion; one batch only"
        ),
        "remaining_open": (
            "584 balanced entries, including 124 split-cubic/split-quintic "
            "entries, plus interval certification and the adaptive lift"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic degree-ten completion-row insertion",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_degree_ten_completion_row_insertion.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 degree-ten completion-row insertion: "
        f"entries={result.closed_entries},"
        f"orbits={result.closed_orbits},"
        f"coefficients={result.left_double_singleton_coefficient:.12g}/"
        f"{result.reversed_double_singleton_coefficient:.12g}/"
        f"{result.double_endpoint_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.degree_ten_inserted.total:.12g},"
        f"margin={result.degree_ten_inserted.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        f"proxy_total={result.remaining_quintic_local_proxy.total:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
