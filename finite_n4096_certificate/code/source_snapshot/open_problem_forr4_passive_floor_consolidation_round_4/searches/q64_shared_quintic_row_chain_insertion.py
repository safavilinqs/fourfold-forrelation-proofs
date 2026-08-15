#!/usr/bin/env python3
"""Shared row/chain contraction for the last 48 q64 quintic entries."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from functools import lru_cache
from json import dumps
from math import comb, prod, sqrt
from pathlib import Path
import sys
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent
    / "open_problem_forr4_passive_floor_consolidation_round_3"
    / "searches"
)
sys.path.insert(0, str(ROUND3_SEARCHES))

from middle_cubic_quintic_pair_contraction import middle_link_maxima  # noqa: E402
import occupation_compatible_sector_optimization as occupation  # noqa: E402
from q64_adjacent_double_cubic_quintic_endpoint_insertion import (  # noqa: E402
    inserted_coefficients as previous_inserted_coefficients,
    remaining_quintic_entries as previous_remaining_quintic_entries,
)
from q64_degree_ten_completion_row_insertion import orbit  # noqa: E402
from q64_paper_target_gate import (  # noqa: E402
    DIMENSION,
    MODES,
    ORDER,
    RESERVE_TARGET,
    THRESHOLD,
    OptimizedLedger,
    ProfileSplit,
    optimize,
)
from q64_post_universal_quintic_gate import quintic_split_depth  # noqa: E402

# Four orientations have a favorable singleton--quintic endpoint followed
# by the split cubic. Three have the same endpoint followed by a whole cubic.
FAVORABLE_ADJACENT_SPLIT = (
    ((1, 5, 3, 3), (0, 1, 2, 3)),
    ((1, 5, 3, 3), (0, 2, 1, 3)),
    ((3, 1, 5, 3), (0, 1, 3, 2)),
    ((3, 1, 5, 3), (0, 1, 4, 1)),
)
FAVORABLE_ADJACENT_WHOLE = (
    ((1, 5, 3, 3), (0, 1, 3, 2)),
    ((1, 5, 3, 3), (0, 2, 3, 1)),
    ((3, 1, 5, 3), (1, 1, 4, 0)),
)

# Here the endpoint singleton lies on the other cubic. A four-record sector
# decomposition retains both cubic--quintic links instead of multiplying
# independent endpoint proxies.
NONFAVORABLE_ADJACENT_SPLIT = (
    ((1, 3, 5, 3), (0, 3, 1, 2)),
    ((1, 3, 5, 3), (0, 3, 2, 1)),
)

# These orientations reduce to a fixed-one singleton--quintic/whole-cubic
# row; the remaining physical link is a unit Schur multiplier.
FIXED_ONE_ADJACENT_WHOLE = (
    ((1, 5, 3, 3), (0, 4, 0, 2)),
    ((3, 1, 5, 3), (1, 1, 1, 3)),
)

# The last orbit has a singleton--whole-cubic endpoint followed by a
# quintic with four fixed cells.
CUBIC_ENDPOINT_FIXED_FOUR = (((3, 1, 3, 5), (1, 1, 0, 4)),)

TARGETS = (
    FAVORABLE_ADJACENT_SPLIT
    + FAVORABLE_ADJACENT_WHOLE
    + NONFAVORABLE_ADJACENT_SPLIT
    + FIXED_ONE_ADJACENT_WHOLE
    + CUBIC_ENDPOINT_FIXED_FOUR
)
RESIDUAL_CLASS_LABELS = (
    "higher_split_only_in_cubic_profile",
    "noncubic_profile",
    "two_split_cubics_one_split_higher",
    "one_split_cubic_no_split_higher",
)
RESIDUAL_CLASS_RESERVE_GATES = (
    0.15571081260134922,
    0.5358557351875237,
    0.34919334312169487,
    0.5571266349301651,
)

IncidenceFamily = Callable[[int], int]


@dataclass(frozen=True)
class SharedQuinticRowChainInsertion:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    closed_orbits: int
    extreme_entries: int
    balanced_entries: int
    fixed_one_record_one_energy: float
    fixed_one_record_three_energy: float
    fixed_one_record_five_energy: float
    fixed_one_total_energy: float
    record_one_middle_maximum: float
    record_three_middle_maximum: float
    favorable_adjacent_split_extreme_coefficient: float
    favorable_adjacent_split_balanced_coefficient: float
    favorable_adjacent_whole_extreme_coefficient: float
    favorable_adjacent_whole_balanced_coefficient: float
    nonfavorable_adjacent_split_extreme_coefficient: float
    nonfavorable_adjacent_split_balanced_coefficient: float
    fixed_one_adjacent_whole_coefficient: float
    cubic_endpoint_fixed_four_coefficient: float
    minimum_coefficient: float
    maximum_coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    remaining_quintic_entries: int
    previous_routing: OptimizedLedger
    inserted_routing: OptimizedLedger
    routing_improvement: float
    adaptive_additive_cap_retaining_allowance: float
    adaptive_multiplier_cap_retaining_allowance: float
    residual_class_labels: tuple[str, ...]
    residual_class_counts: tuple[int, ...]
    residual_class_frozen_targets: tuple[float, ...]
    residual_class_reserve_gates: tuple[float, ...]


def shared_quintic_entries() -> tuple[ProfileSplit, ...]:
    """Return the twelve complement/reversal orbits closed here."""

    return tuple(sorted({entry for target in TARGETS for entry in orbit(target)}))


def endpoint_quintic_fixed_one_record_energies(
    order: int,
) -> tuple[float, float, float]:
    """Split the exact fixed-one endpoint energy by the other record.

    The singleton link forces one quintic record to have size one. The
    return value gives the squared endpoint-moment sums for other record
    sizes one, three, and five. The formulas refine the exact ``F_1`` slice
    by the column patterns ``5``, ``4+1``, ``3+2``, and ``2+2+1``.
    """

    q = order
    if q < 4 or q & (q - 1):
        raise ValueError(("power-of-two order at least four required", q))
    w0 = 1 / q**2
    w1 = 1 / (q**2 * (q - 1) ** 2)
    w2 = 4 / (q**2 * (q - 1) ** 2 * (q - 2) ** 2)
    totals = {1: 0.0, 3: 0.0, 5: 0.0}

    # Column pattern 5: five distinct rows, hence other record five.
    totals[5] += q * comb(q, 5) * w0

    # Column pattern 4+1. The singleton row either belongs to the four-set
    # (record three) or does not (record five). Exactly Z_4 four-sets have
    # zero xor and receive w0 rather than w1.
    column_choices = q * (q - 1)
    zero_xor_four_sets = q * (q - 1) * (q - 2) // 24
    four_sets = comb(q, 4)
    totals[3] += column_choices * (
        4 * zero_xor_four_sets * w0
        + 4 * (four_sets - zero_xor_four_sets) * w1
    )
    totals[5] += column_choices * (
        (q - 4) * zero_xor_four_sets * w0
        + (q - 4) * (four_sets - zero_xor_four_sets) * w1
    )

    # Column pattern 3+2. If the row sets overlap in t cells, the other
    # parity record has size 5-2t.
    for overlap in range(3):
        count = (
            q
            * (q - 1)
            * comb(q, 3)
            * comb(3, overlap)
            * comb(q - 3, 2 - overlap)
        )
        totals[5 - 2 * overlap] += count * w1

    # Column pattern 2+2+1. Classify the two row pairs by symmetric-
    # difference size 0, 2, or 4. Equal pair xors receive w1; all other
    # pair choices receive w2.
    column_choices = q * comb(q - 1, 2)
    pairs = comb(q, 2)
    total_counts = {
        1: pairs * q + pairs * 2 * (q - 2) * 2,
        3: pairs * 2 * (q - 2) * (q - 2)
        + pairs * comb(q - 2, 2) * 4,
        5: pairs * comb(q - 2, 2) * (q - 4),
    }
    equal_xor_counts = {
        1: pairs * q,
        3: pairs * (q // 2 - 1) * 4,
        5: pairs * (q // 2 - 1) * (q - 4),
    }
    for record in (1, 3, 5):
        totals[record] += column_choices * (
            equal_xor_counts[record] * w1
            + (total_counts[record] - equal_xor_counts[record]) * w2
        )

    # Cell transitivity converts each global sector total to its fixed-one
    # slice: every support has five cells among q^2 possible fixed cells.
    return tuple(5 * totals[record] / q**2 for record in (1, 3, 5))


def integer_partitions(total: int, maximum: int | None = None):
    """Yield decreasing positive partitions of a small integer."""

    if total == 0:
        yield ()
        return
    if maximum is None or maximum > total:
        maximum = total
    for first in range(maximum, 0, -1):
        for tail in integer_partitions(total - first, first):
            yield (first,) + tail


def odd_record_extension_count(
    order: int,
    fixed_counts: tuple[int, ...],
    degree: int,
    record: int,
) -> int:
    """Count one-axis record-``record`` extensions of a fixed pattern."""

    fixed_size = sum(fixed_counts)
    groups = fixed_counts + (0,) * (order - len(fixed_counts))
    dynamic = {(0, 0): 1}
    for fixed in groups:
        updated: dict[tuple[int, int], int] = {}
        for (added, odd), count in dynamic.items():
            remaining = degree - fixed_size - added
            for extra in range(min(order - fixed, remaining) + 1):
                key = (added + extra, odd + (fixed + extra) % 2)
                updated[key] = updated.get(key, 0) + count * comb(
                    order - fixed, extra
                )
        dynamic = updated
    return dynamic.get((degree - fixed_size, record), 0)


@lru_cache(maxsize=None)
def odd_record_incidence(
    order: int,
    degree: int,
    record: int,
    selected: int,
) -> int:
    """Maximum one-axis incidence with a prescribed odd-record size."""

    return max(
        odd_record_extension_count(order, partition, degree, record)
        for partition in integer_partitions(selected)
        if len(partition) <= order
    )


def two_axis_relaxed_incidence(
    order: int,
    degree: int,
    left_record: int,
    right_record: int,
    selected: int,
) -> int:
    """Safe incidence for a bidegree sector by intersecting relaxations."""

    return min(
        odd_record_incidence(order, degree, left_record, selected),
        odd_record_incidence(order, degree, right_record, selected),
    )


def rank_incidence_coefficient(
    profile: tuple[int, ...],
    split: tuple[int, ...],
    families: tuple[IncidenceFamily, ...],
    maximum_entry: float,
    order: int,
) -> float:
    """Combine the cut-rank and row/column incidence bounds."""

    complement = tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )
    row = prod(
        family(selected)
        for family, selected in zip(families, split, strict=True)
    )
    column = prod(
        family(selected)
        for family, selected in zip(families, complement, strict=True)
    )
    rank = order ** min(sum(split), sum(complement)) * maximum_entry
    return min(rank, maximum_entry * sqrt(row), maximum_entry * sqrt(column))


def favorable_adjacent_split_coefficients(
    order: int = ORDER,
) -> tuple[float, float]:
    """Record-compatible complete-row bounds for the first four orbits."""

    q = order
    dimension = q * q
    record_one, record_three, universal = middle_link_maxima(q)
    fixed_one = endpoint_quintic_fixed_one_record_energies(q)
    quintic = occupation.endpoint_quintic_singleton_slice_energies(q)
    extreme = sqrt(
        dimension
        * (
            fixed_one[0] * record_one**2
            + fixed_one[1] * record_three**2
        )
    )
    # For the balanced cut, the exact fixed-three endpoint slice is already
    # comfortably below the gate; the universal record-one maximum is safe.
    balanced = sqrt(dimension * quintic[3] * universal**2)
    return extreme, balanced


def favorable_adjacent_whole_coefficients(
    order: int = ORDER,
) -> tuple[float, float]:
    """Endpoint/whole-cubic complete-row bounds for three orbits."""

    _, _, universal = middle_link_maxima(order)
    quintic = occupation.endpoint_quintic_singleton_slice_energies(order)
    return sqrt(quintic[4]) * universal, sqrt(quintic[3]) * universal


def nonfavorable_adjacent_split_coefficient(
    target: ProfileSplit,
    order: int = ORDER,
) -> float:
    """Four-sector incidence bound for the endpoint-cubic/quintic chain."""

    profile, split = target
    if profile != (1, 3, 5, 3):
        raise ValueError(("canonical profile required", target))
    q = order
    maxima = dict(zip((1, 3), middle_link_maxima(q)[:2], strict=True))

    def singleton(selected: int) -> int:
        return q * q if selected == 0 else 1

    coefficient = 0.0
    for first_record in (1, 3):
        for second_record in (1, 3):
            endpoint = (
                1 / (q * (q - 1))
                if first_record == 1
                else 1 / q
            )
            maximum_entry = (
                endpoint * maxima[first_record] * maxima[second_record]
            )
            families: tuple[IncidenceFamily, ...] = (
                singleton,
                lambda selected, record=first_record: (
                    two_axis_relaxed_incidence(
                        q, 3, 1, record, selected
                    )
                ),
                lambda selected, left=first_record, right=second_record: (
                    two_axis_relaxed_incidence(
                        q, 5, left, right, selected
                    )
                ),
                lambda selected, record=second_record: odd_record_incidence(
                    q, 3, record, selected
                ),
            )
            coefficient += rank_incidence_coefficient(
                profile, split, families, maximum_entry, q
            )
    return coefficient


def fixed_one_adjacent_whole_coefficient(order: int = ORDER) -> float:
    """Fixed-one endpoint energy followed by one whole cubic."""

    _, _, universal = middle_link_maxima(order)
    fixed_one = occupation.endpoint_quintic_singleton_slice_energies(order)[1]
    return sqrt(fixed_one) * universal


def cubic_endpoint_fixed_four_coefficient(order: int = ORDER) -> float:
    """Singleton--whole-cubic row followed by a fixed-four quintic."""

    q = order
    record_one, record_three, _ = middle_link_maxima(q)
    record_one_term = (
        record_one
        * sqrt(odd_record_incidence(q, 5, 1, 4))
        / (q * (q - 1))
    )
    record_three_term = (
        record_three
        * sqrt(odd_record_incidence(q, 5, 3, 4))
        / q
    )
    return record_one_term + record_three_term


def coefficient_map(order: int = ORDER) -> dict[ProfileSplit, float]:
    """Return the proved coefficient for every closed entry."""

    adjacent_extreme, adjacent_balanced = (
        favorable_adjacent_split_coefficients(order)
    )
    whole_extreme, whole_balanced = favorable_adjacent_whole_coefficients(
        order
    )
    fixed_one_whole = fixed_one_adjacent_whole_coefficient(order)
    cubic_fixed_four = cubic_endpoint_fixed_four_coefficient(order)
    result: dict[ProfileSplit, float] = {}

    def insert(target: ProfileSplit, value: float) -> None:
        for entry in orbit(target):
            if entry in result:
                raise AssertionError(("shared-quintic orbit overlap", entry))
            result[entry] = value

    for target in FAVORABLE_ADJACENT_SPLIT:
        insert(
            target,
            adjacent_extreme
            if quintic_split_depth(target) == 1
            else adjacent_balanced,
        )
    for target in FAVORABLE_ADJACENT_WHOLE:
        insert(
            target,
            whole_extreme
            if quintic_split_depth(target) == 1
            else whole_balanced,
        )
    for target in NONFAVORABLE_ADJACENT_SPLIT:
        insert(target, nonfavorable_adjacent_split_coefficient(target, order))
    for target in FIXED_ONE_ADJACENT_WHOLE:
        insert(target, fixed_one_whole)
    for target in CUBIC_ENDPOINT_FIXED_FOUR:
        insert(target, cubic_fixed_four)
    return result


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = previous_inserted_coefficients()
    result.update(coefficient_map())
    return result


def remaining_quintic_entries() -> tuple[ProfileSplit, ...]:
    closed = set(shared_quintic_entries())
    return tuple(
        entry
        for entry in previous_remaining_quintic_entries()
        if entry not in closed
    )


def diagnostic() -> SharedQuinticRowChainInsertion:
    from q64_remaining_class_gates import partition_remaining

    entries = shared_quintic_entries()
    fixed_one = endpoint_quintic_fixed_one_record_energies(ORDER)
    record_one, record_three, _ = middle_link_maxima(ORDER)
    adjacent_extreme, adjacent_balanced = (
        favorable_adjacent_split_coefficients()
    )
    whole_extreme, whole_balanced = favorable_adjacent_whole_coefficients()
    nonfavorable = tuple(
        nonfavorable_adjacent_split_coefficient(target)
        for target in NONFAVORABLE_ADJACENT_SPLIT
    )
    fixed_one_whole = fixed_one_adjacent_whole_coefficient()
    cubic_fixed_four = cubic_endpoint_fixed_four_coefficient()
    coefficients = coefficient_map()
    previous = optimize(mapped_coefficients=previous_inserted_coefficients())
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    previous_proved = 380
    partition = partition_remaining()
    residual_entries = tuple(
        partition[label] for label in RESIDUAL_CLASS_LABELS
    )
    inserted_coefficients_map = inserted_coefficients()
    return SharedQuinticRowChainInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        closed_orbits=len(TARGETS),
        extreme_entries=sum(quintic_split_depth(entry) == 1 for entry in entries),
        balanced_entries=sum(
            quintic_split_depth(entry) == 2 for entry in entries
        ),
        fixed_one_record_one_energy=fixed_one[0],
        fixed_one_record_three_energy=fixed_one[1],
        fixed_one_record_five_energy=fixed_one[2],
        fixed_one_total_energy=sum(fixed_one),
        record_one_middle_maximum=record_one,
        record_three_middle_maximum=record_three,
        favorable_adjacent_split_extreme_coefficient=adjacent_extreme,
        favorable_adjacent_split_balanced_coefficient=adjacent_balanced,
        favorable_adjacent_whole_extreme_coefficient=whole_extreme,
        favorable_adjacent_whole_balanced_coefficient=whole_balanced,
        nonfavorable_adjacent_split_extreme_coefficient=nonfavorable[0],
        nonfavorable_adjacent_split_balanced_coefficient=nonfavorable[1],
        fixed_one_adjacent_whole_coefficient=fixed_one_whole,
        cubic_endpoint_fixed_four_coefficient=cubic_fixed_four,
        minimum_coefficient=min(coefficients.values()),
        maximum_coefficient=max(coefficients.values()),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved + len(entries),
        remaining_open_entries=888 - previous_proved - len(entries),
        remaining_quintic_entries=len(remaining_quintic_entries()),
        previous_routing=previous,
        inserted_routing=inserted,
        routing_improvement=previous.total - inserted.total,
        adaptive_additive_cap_retaining_allowance=(
            inserted.margin_to_one_third - RESERVE_TARGET
        ),
        adaptive_multiplier_cap_retaining_allowance=(
            (THRESHOLD - RESERVE_TARGET) / inserted.total
        ),
        residual_class_labels=RESIDUAL_CLASS_LABELS,
        residual_class_counts=tuple(len(entries) for entries in residual_entries),
        residual_class_frozen_targets=tuple(
            inserted_coefficients_map[entries[0]]
            for entries in residual_entries
        ),
        residual_class_reserve_gates=RESIDUAL_CLASS_RESERVE_GATES,
    )


def artifact_text(result: SharedQuinticRowChainInsertion) -> str:
    payload = {
        "schema": "round4_q64_shared_quintic_row_chain_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal row/chain theorem for all 48 remaining q64 "
            "split-cubic/split-quintic entries; exact record-compatible "
            "fixed-one endpoint decomposition, finite sector list, and "
            "rank/incidence bounds; floating q64 Perron insertion; one "
            "batch only"
        ),
        "adaptive_requirement": (
            "the displayed additive and multiplicative caps are algebraic "
            "requirements conditional on closing the other 460 entries at "
            "their frozen coefficient targets; they are not an adaptive "
            "recurrence or passive theorem"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic shared-quintic insertion artifact",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_shared_quintic_row_chain_insertion.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 shared quintic row/chain insertion: "
        f"entries={result.closed_entries},"
        f"orbits={result.closed_orbits},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.inserted_routing.total:.12g},"
        f"margin={result.inserted_routing.margin_to_one_third:.12g},"
        "adaptive_additive_cap="
        f"{result.adaptive_additive_cap_retaining_allowance:.12g},"
        f"remaining_open={result.remaining_open_entries},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
