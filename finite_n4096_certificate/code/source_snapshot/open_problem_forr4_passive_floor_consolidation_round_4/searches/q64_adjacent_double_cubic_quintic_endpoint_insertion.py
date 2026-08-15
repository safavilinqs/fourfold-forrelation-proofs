#!/usr/bin/env python3
"""Adjacent double-cubic plus quintic-endpoint insertion at q64."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from math import comb, sqrt
from pathlib import Path

from q64_balanced_pair_triple_mask_insertion import (
    inserted_coefficients as pair_triple_inserted_coefficients,
    pair_triple_disjointness_factor,
    remaining_quintic_entries as pre_adjacent_double_cubic_entries,
)
from q64_degree_ten_completion_row_insertion import orbit
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
from q64_dual_endpoint_schur_insertion import has_favorable_quintic_singleton
from leading_balanced_disjointness_contraction import (
    disjointness_schur_factor,
)
import occupation_compatible_sector_optimization as occupation


ROOT = Path(__file__).resolve().parents[1]
WHOLE_ENDPOINT_EXTREME = ((3, 3, 1, 5), (0, 1, 1, 4))
WHOLE_MIDDLE_EXTREME = ((3, 3, 1, 5), (1, 0, 1, 4))
WHOLE_MIDDLE_BALANCED = ((3, 3, 1, 5), (1, 3, 0, 2))
GENERIC_S_C_W_Q_EXTREME = ((1, 3, 3, 5), (0, 2, 0, 4))
GENERIC_S_C_W_Q_COMPLEMENT_EXTREME = ((1, 3, 3, 5), (0, 2, 3, 1))
GENERIC_S_W_C_Q_BALANCED = ((1, 3, 3, 5), (0, 3, 1, 2))
GENERIC_S_W_C_Q_EXTREME = ((1, 3, 3, 5), (0, 3, 2, 1))
GENERIC_W_C_S_Q_EXTREME = ((3, 3, 1, 5), (0, 2, 0, 4))
TARGETS = (
    WHOLE_ENDPOINT_EXTREME,
    WHOLE_MIDDLE_EXTREME,
    WHOLE_MIDDLE_BALANCED,
    GENERIC_S_C_W_Q_EXTREME,
    GENERIC_S_C_W_Q_COMPLEMENT_EXTREME,
    GENERIC_S_W_C_Q_BALANCED,
    GENERIC_S_W_C_Q_EXTREME,
    GENERIC_W_C_S_Q_EXTREME,
)


@dataclass(frozen=True)
class AdjacentDoubleCubicQuinticEndpointInsertion:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    closed_orbits: int
    extreme_entries: int
    balanced_entries: int
    endpoint_middle_split_pairs: tuple[tuple[int, int], ...]
    forward_double_cubic_coefficient: float
    complement_double_cubic_coefficient: float
    double_cubic_coefficient: float
    quintic_fixed_four_energy: float
    quintic_endpoint_factor: float
    coefficient: float
    middle_whole_extreme_coefficient: float
    quintic_fixed_three_energy: float
    middle_whole_balanced_coefficient: float
    generic_extreme_mask_factor: float
    generic_balanced_mask_factor: float
    minimum_coefficient: float
    maximum_coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    inserted_routing: OptimizedLedger
    routing_change: float
    remaining_quintic_entries: int
    remaining_extreme_entries: int
    remaining_balanced_entries: int
    remaining_quintic_local_proxy: OptimizedLedger
    proxy_reserve_after_declared_allowance: float


def adjacent_double_cubic_entries() -> tuple[ProfileSplit, ...]:
    return tuple(sorted({entry for target in TARGETS for entry in orbit(target)}))


def endpoint_middle_split_pair(entry: ProfileSplit) -> tuple[int, int]:
    profile, split = entry
    if profile[:3] == (3, 3, 1):
        return split[0], split[1]
    if profile[1:] == (1, 3, 3):
        return split[3], split[2]
    if profile[:3] == (1, 3, 3):
        return split[2], split[1]
    if profile[1:] == (3, 3, 1):
        return split[1], split[2]
    raise ValueError(("missing whole-cubic/split-cubic/singleton chain", entry))


def uniform_double_cubic_coefficient(
    endpoint_size: int,
    middle_size: int,
    order: int = ORDER,
) -> float:
    """Inherited arbitrary-law incidence bound for ``M33 M31``."""

    maximum_squared_entry = (order + 2) ** 2 / (
        order * order * (order - 1) ** 2 * (order - 2) ** 2
    )
    endpoint_degrees = (
        order * comb(order, 3)
        + order * order * (order - 1) * comb(order, 2),
        comb(order - 1, 2)
        + (order - 1) * comb(order, 2)
        + order * (order - 1) ** 2,
        order * order - 2,
        1,
    )
    middle_degrees = (
        order * order * (order - 1) ** 2,
        3 * (order - 1) ** 2,
        2 * (order - 1),
        1,
    )
    row = sqrt(
        maximum_squared_entry
        * endpoint_degrees[endpoint_size]
        * middle_degrees[middle_size]
    ) / (order - 1)
    column = sqrt(
        maximum_squared_entry
        * endpoint_degrees[3 - endpoint_size]
        * middle_degrees[3 - middle_size]
    ) / (order * (order - 1))
    return min(row, column)


def coefficient(order: int = ORDER) -> float:
    double_cubic = max(
        uniform_double_cubic_coefficient(0, 1, order),
        uniform_double_cubic_coefficient(3, 2, order),
    )
    quintic = occupation.endpoint_quintic_singleton_slice_energies(order)[4]
    return double_cubic * sqrt(quintic)


def target_coefficient(target: ProfileSplit, order: int = ORDER) -> float:
    pairs = {endpoint_middle_split_pair(entry) for entry in orbit(target)}
    double_cubic = max(
        uniform_double_cubic_coefficient(*pair, order=order) for pair in pairs
    )
    quintic_depth = quintic_split_depth(target)
    if has_favorable_quintic_singleton(target):
        fixed_count = 4 if quintic_depth == 1 else 3
        quintic = occupation.endpoint_quintic_singleton_slice_energies(order)[
            fixed_count
        ]
        quintic_factor = sqrt(quintic)
    elif quintic_depth == 1:
        quintic_factor = disjointness_schur_factor(order * order, 4)
    else:
        quintic_factor = pair_triple_disjointness_factor()
    return double_cubic * quintic_factor


def coefficient_map(order: int = ORDER) -> dict[ProfileSplit, float]:
    result: dict[ProfileSplit, float] = {}
    for target in TARGETS:
        value = target_coefficient(target, order)
        for entry in orbit(target):
            if entry in result:
                raise AssertionError(("adjacent-double-cubic overlap", entry))
            result[entry] = value
    return result


def remaining_quintic_entries() -> tuple[ProfileSplit, ...]:
    closed = set(adjacent_double_cubic_entries())
    return tuple(
        entry
        for entry in pre_adjacent_double_cubic_entries()
        if entry not in closed
    )


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = pair_triple_inserted_coefficients()
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


def diagnostic() -> AdjacentDoubleCubicQuinticEndpointInsertion:
    entries = adjacent_double_cubic_entries()
    remaining = remaining_quintic_entries()
    forward = uniform_double_cubic_coefficient(0, 1)
    complement = uniform_double_cubic_coefficient(3, 2)
    quintic_energy = occupation.endpoint_quintic_singleton_slice_energies(
        ORDER
    )[4]
    quintic_three = occupation.endpoint_quintic_singleton_slice_energies(
        ORDER
    )[3]
    previous = optimize(
        mapped_coefficients=pair_triple_inserted_coefficients()
    )
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    proxy = optimize(
        mapped_coefficients=remaining_quintic_local_proxy_coefficients()
    )
    previous_proved = 348
    coefficients = coefficient_map()
    return AdjacentDoubleCubicQuinticEndpointInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        closed_orbits=len({frozenset(orbit(entry)) for entry in entries}),
        extreme_entries=sum(quintic_split_depth(entry) == 1 for entry in entries),
        balanced_entries=sum(
            quintic_split_depth(entry) == 2 for entry in entries
        ),
        endpoint_middle_split_pairs=tuple(
            sorted(endpoint_middle_split_pair(entry) for entry in entries)
        ),
        forward_double_cubic_coefficient=forward,
        complement_double_cubic_coefficient=complement,
        double_cubic_coefficient=max(forward, complement),
        quintic_fixed_four_energy=quintic_energy,
        quintic_endpoint_factor=sqrt(quintic_energy),
        coefficient=coefficient(),
        middle_whole_extreme_coefficient=target_coefficient(
            WHOLE_MIDDLE_EXTREME
        ),
        quintic_fixed_three_energy=quintic_three,
        middle_whole_balanced_coefficient=target_coefficient(
            WHOLE_MIDDLE_BALANCED
        ),
        generic_extreme_mask_factor=disjointness_schur_factor(DIMENSION, 4),
        generic_balanced_mask_factor=pair_triple_disjointness_factor(),
        minimum_coefficient=min(coefficients.values()),
        maximum_coefficient=max(coefficients.values()),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        inserted_routing=inserted,
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


def artifact_text(
    result: AdjacentDoubleCubicQuinticEndpointInsertion,
) -> str:
    payload = {
        "schema": (
            "round4_q64_adjacent_double_cubic_quintic_endpoint_insertion_v1"
        ),
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal adjacent M33-M31 incidence contraction "
            "composed with the exact favorable fixed-four quintic endpoint "
            "Schur factor; floating q64 Perron insertion; one batch only"
        ),
        "remaining_open": (
            "508 balanced entries, including 48 degree-twelve "
            "split-cubic/split-quintic entries, plus interval certification "
            "and the adaptive lift"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic adjacent-double-cubic insertion",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = (
            ROOT
            / "artifacts"
            / "q64_adjacent_double_cubic_quintic_endpoint_insertion.json"
        )
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 adjacent double-cubic quintic-endpoint insertion: "
        f"entries={result.closed_entries},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.inserted_routing.total:.12g},"
        f"margin={result.inserted_routing.margin_to_one_third:.12g},"
        f"remaining_quintic={result.remaining_quintic_entries},"
        f"proxy_total={result.remaining_quintic_local_proxy.total:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
