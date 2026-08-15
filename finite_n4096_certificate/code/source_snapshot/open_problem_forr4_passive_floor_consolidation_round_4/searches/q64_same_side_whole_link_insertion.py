#!/usr/bin/env python3
"""Close 96 q64 entries through one same-side whole link and one mask."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from json import dumps
from math import comb, sqrt
from pathlib import Path

from q64_degree_ten_completion_row_insertion import orbit
from q64_paper_target_gate import (
    DIMENSION,
    MODES,
    ORDER,
    RESERVE_TARGET,
    THRESHOLD,
    OptimizedLedger,
    ProfileSplit,
    optimize,
)
from q64_remaining_class_gates import partition_remaining
from q64_whole_higher_split_cubic_insertion import (
    inserted_coefficients as previous_inserted_coefficients,
)


ROOT = Path(__file__).resolve().parents[1]
CLASS_LABEL = "higher_split_only_in_cubic_profile"
REMAINING_CLASS_RESERVE_GATE = 0.190775718804


@dataclass(frozen=True)
class SameSideWholeLinkInsertion:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    closed_orbits: int
    singleton_singleton_entries: int
    singleton_cubic_entries: int
    degree_five_extreme_entries: int
    degree_five_balanced_entries: int
    degree_seven_extreme_entries: int
    degree_seven_two_five_entries: int
    degree_seven_three_four_entries: int
    degree_five_extreme_mask_factor: float
    degree_five_balanced_mask_factor: float
    degree_seven_extreme_mask_factor: float
    degree_seven_two_five_mask_factor: float
    degree_seven_three_four_mask_factor: float
    link_factor: float
    minimum_coefficient: float
    maximum_coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    inserted_routing: OptimizedLedger
    routing_improvement: float
    reserve_after_declared_allowance: float
    adaptive_multiplier_cap_retaining_allowance: float
    remaining_class_reserve_gate: float


def higher_split(entry: ProfileSplit) -> tuple[int, int, int]:
    """Return the index, degree, and smaller side of the unique high split."""

    profile, split = entry
    values = tuple(
        (index, degree, min(selected, degree - selected))
        for index, (degree, selected) in enumerate(
            zip(profile, split, strict=True)
        )
        if degree > 3 and 0 < selected < degree
    )
    if len(values) != 1:
        raise ValueError(("unique split higher block required", entry))
    return values[0]


def same_side_whole_links(entry: ProfileSplit) -> tuple[tuple[int, int], ...]:
    """Adjacent whole-block degree pairs lying on the same cut side."""

    profile, split = entry
    result = []
    for index in range(3):
        left_whole = split[index] in (0, profile[index])
        right_whole = split[index + 1] in (0, profile[index + 1])
        same_side = (split[index] == profile[index]) == (
            split[index + 1] == profile[index + 1]
        )
        if left_whole and right_whole and same_side:
            result.append((profile[index], profile[index + 1]))
    return tuple(result)


def subset_disjointness_factor(left: int, right: int) -> float:
    """Direct-sum Schur factor for disjoint left- and right-subsets."""

    if left < 0 or right < 0:
        raise ValueError(("nonnegative subset sizes required", left, right))
    return sum(
        sqrt(comb(left, level) * comb(right, level))
        for level in range(min(left, right) + 1)
    )


def target_entries() -> tuple[ProfileSplit, ...]:
    return tuple(
        entry
        for entry in partition_remaining()[CLASS_LABEL]
        if same_side_whole_links(entry)
    )


def remaining_entries() -> tuple[ProfileSplit, ...]:
    closed = set(target_entries())
    return tuple(
        entry
        for entry in partition_remaining()[CLASS_LABEL]
        if entry not in closed
    )


def entry_coefficient(entry: ProfileSplit, order: int = ORDER) -> float:
    """One whole-link factor times the unique split-block mask factor."""

    _, degree, depth = higher_split(entry)
    return subset_disjointness_factor(depth, degree - depth) / order


def coefficient_map() -> dict[ProfileSplit, float]:
    return {entry: entry_coefficient(entry) for entry in target_entries()}


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = previous_inserted_coefficients()
    result.update(coefficient_map())
    return result


def diagnostic() -> SameSideWholeLinkInsertion:
    entries = target_entries()
    coefficients = coefficient_map()
    shapes = Counter(higher_split(entry)[1:] for entry in entries)
    link_types = Counter(
        "singleton_singleton"
        if (1, 1) in same_side_whole_links(entry)
        else "singleton_cubic"
        for entry in entries
    )
    previous = optimize(mapped_coefficients=previous_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    previous_proved = 712
    return SameSideWholeLinkInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        closed_orbits=len({frozenset(orbit(entry)) for entry in entries}),
        singleton_singleton_entries=link_types["singleton_singleton"],
        singleton_cubic_entries=link_types["singleton_cubic"],
        degree_five_extreme_entries=shapes[(5, 1)],
        degree_five_balanced_entries=shapes[(5, 2)],
        degree_seven_extreme_entries=shapes[(7, 1)],
        degree_seven_two_five_entries=shapes[(7, 2)],
        degree_seven_three_four_entries=shapes[(7, 3)],
        degree_five_extreme_mask_factor=subset_disjointness_factor(1, 4),
        degree_five_balanced_mask_factor=subset_disjointness_factor(2, 3),
        degree_seven_extreme_mask_factor=subset_disjointness_factor(1, 6),
        degree_seven_two_five_mask_factor=subset_disjointness_factor(2, 5),
        degree_seven_three_four_mask_factor=subset_disjointness_factor(3, 4),
        link_factor=1 / ORDER,
        minimum_coefficient=min(coefficients.values()),
        maximum_coefficient=max(coefficients.values()),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        previous_routing=previous,
        inserted_routing=inserted,
        routing_improvement=previous.total - inserted.total,
        reserve_after_declared_allowance=(
            inserted.margin_to_one_third - RESERVE_TARGET
        ),
        adaptive_multiplier_cap_retaining_allowance=(
            (THRESHOLD - RESERVE_TARGET) / inserted.total
        ),
        remaining_class_reserve_gate=REMAINING_CLASS_RESERVE_GATE,
    )


def artifact_text(result: SameSideWholeLinkInsertion) -> str:
    payload = {
        "schema": "round4_q64_same_side_whole_link_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal q64 theorem for 96 entries: a same-side "
            "whole singleton-singleton or singleton-cubic link supplies "
            "1/q, the only split-block disjointness mask is restored by "
            "an exact inclusion direct-sum Schur factor, and all completed "
            "physical links are unit cross-Gram multipliers; floating "
            "Perron insertion; one batch only"
        ),
        "adaptive_requirement": (
            "the local 96-entry theorem survives; its cumulative count, "
            "routing reserve, and adaptive allowance are withdrawn because "
            "the upstream universal insertion is quarantined"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic same-side whole-link artifact",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_same_side_whole_link_insertion.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 same-side whole-link insertion: "
        f"entries={result.closed_entries},"
        f"links={result.singleton_singleton_entries}/"
        f"{result.singleton_cubic_entries},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.inserted_routing.total:.12g},"
        f"margin={result.inserted_routing.margin_to_one_third:.12g},"
        f"remaining_open={result.remaining_open_entries},"
        f"remaining_gate={result.remaining_class_reserve_gate:.12g},"
        "status=local_96_entry_theorem_cumulative_values_withdrawn"
    )


if __name__ == "__main__":
    main()
