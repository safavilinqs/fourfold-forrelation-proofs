#!/usr/bin/env python3
"""Complete-row contraction for the q64 adjacent fixed-pair orbit.

For profile (1,1,3,5) and split (0,1,2,2), a row fixes the second
singleton, two cubic cells, and two quintic cells.  Keeping M_13 M_35 as one
Schur feature exposes the fact that only O(q) cubic completions can carry
the largest endpoint amplitudes.
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
    RESERVE_TARGET,
    OptimizedLedger,
    ProfileSplit,
    optimize,
)
import occupation_compatible_sector_optimization as occupation
from q64_fixed_singleton_pair_contraction import (
    fixed_singleton_pair_entries,
)
from q64_post_universal_quintic_gate import (
    quintic_entries,
    quintic_split_depth,
)
from q64_reversed_middle_pair_contraction import reversed_middle_pair_entries
from q64_shifted_middle_pair_contraction import shifted_middle_pair_entries
from q64_universal_double_cubic_insertion import (
    inserted_coefficients as double_cubic_inserted_coefficients,
)
from whole_cubic_middle_pair_contraction import record_sector_bounds


ROOT = Path(__file__).resolve().parents[1]
PROFILE = (1, 1, 3, 5)
SPLIT = (0, 1, 2, 2)


@dataclass(frozen=True)
class FixedPairAdjacentRowContraction:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    record_one_cubic_completions: int
    record_three_high_cubic_completions: int
    record_three_low_cubic_completions: int
    record_one_middle_squared: float
    record_three_high_middle_squared: float
    record_three_low_middle_squared: float
    record_one_tail_bound: float
    record_three_same_pair_tail_bound: float
    record_three_distinct_pair_tail_bound: float
    record_three_tail_bound: float
    record_one_row_energy_bound: float
    record_three_row_energy_bound: float
    total_row_energy_bound: float
    coefficient: float
    previous_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    fixed_pair_inserted: OptimizedLedger
    margin_improvement: float
    remaining_quintic_proxy_entries: int
    remaining_quintic_local_proxy: OptimizedLedger
    proxy_reserve_after_declared_allowance: float


def fixed_pair_adjacent_entries() -> tuple[ProfileSplit, ...]:
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


def cubic_completion_counts(order: int = ORDER) -> tuple[int, int, int]:
    """Maximum low/high completion counts through a fixed cubic pair.

    The classes are an L-shaped record-one cubic, a horizontal record-three
    cubic with endpoint amplitude 1/q, and every remaining record-three
    cubic with endpoint amplitude 1/[q(q-1)].
    """

    if order < 4 or order & (order - 1):
        raise ValueError(("power-of-two order at least four required", order))
    return 2 * (order - 1), order - 2, (order - 1) * (order - 2)


def fixed_pair_adjacent_energy_parts(
    order: int = ORDER,
) -> tuple[float, float]:
    """Return record-one and record-three complete-row energy bounds."""

    q = order
    record_one_count, record_three_high, record_three_low = (
        cubic_completion_counts(q)
    )
    bounds = record_sector_bounds(q)
    record_one_tail = bounds[2]
    record_three_tail = max(bounds[5], bounds[6])
    low_middle_squared = 1 / (q * q * (q - 1) ** 2)
    record_one = record_one_count * low_middle_squared * record_one_tail
    record_three = record_three_tail * (
        record_three_high / q**2
        + record_three_low * low_middle_squared
    )
    return record_one, record_three


def fixed_pair_adjacent_coefficient(order: int = ORDER) -> float:
    return sqrt(sum(fixed_pair_adjacent_energy_parts(order)))


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = double_cubic_inserted_coefficients()
    coefficient = fixed_pair_adjacent_coefficient()
    for entry in fixed_pair_adjacent_entries():
        result[entry] = coefficient
    return result


def remaining_quintic_local_proxy_coefficients() -> dict[ProfileSplit, float]:
    """Assign inherited local scales only to still-open quintic entries."""

    result = inserted_coefficients()
    cubic_slice = occupation.endpoint_singleton_slice_energies(ORDER)[2]
    quintic_slices = occupation.endpoint_quintic_singleton_slice_energies(
        ORDER
    )
    extreme = sqrt(cubic_slice * quintic_slices[4])
    balanced = sqrt(cubic_slice * quintic_slices[3])
    closed = (
        set(shifted_middle_pair_entries())
        | set(reversed_middle_pair_entries())
        | set(fixed_singleton_pair_entries())
        | set(fixed_pair_adjacent_entries())
    )
    for entry in quintic_entries():
        if entry not in closed:
            result[entry] = (
                extreme if quintic_split_depth(entry) == 1 else balanced
            )
    return result


def diagnostic() -> FixedPairAdjacentRowContraction:
    q = ORDER
    counts = cubic_completion_counts(q)
    bounds = record_sector_bounds(q)
    record_one, record_three = fixed_pair_adjacent_energy_parts(q)
    row_energy = record_one + record_three
    previous = optimize(
        mapped_coefficients=double_cubic_inserted_coefficients()
    )
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    local_proxy = optimize(
        mapped_coefficients=remaining_quintic_local_proxy_coefficients()
    )
    previous_proved = 260
    remaining_quintic = len(quintic_entries()) - 16 - len(
        fixed_pair_adjacent_entries()
    )
    return FixedPairAdjacentRowContraction(
        order=q,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(fixed_pair_adjacent_entries()),
        record_one_cubic_completions=counts[0],
        record_three_high_cubic_completions=counts[1],
        record_three_low_cubic_completions=counts[2],
        record_one_middle_squared=1 / (q * q * (q - 1) ** 2),
        record_three_high_middle_squared=1 / q**2,
        record_three_low_middle_squared=1 / (q * q * (q - 1) ** 2),
        record_one_tail_bound=bounds[2],
        record_three_same_pair_tail_bound=bounds[5],
        record_three_distinct_pair_tail_bound=bounds[6],
        record_three_tail_bound=max(bounds[5], bounds[6]),
        record_one_row_energy_bound=record_one,
        record_three_row_energy_bound=record_three,
        total_row_energy_bound=row_energy,
        coefficient=sqrt(row_energy),
        previous_proved_entries=previous_proved,
        total_proved_entries=previous_proved
        + len(fixed_pair_adjacent_entries()),
        remaining_open_entries=(
            888 - previous_proved - len(fixed_pair_adjacent_entries())
        ),
        previous_routing=previous,
        fixed_pair_inserted=inserted,
        margin_improvement=(
            inserted.margin_to_one_third - previous.margin_to_one_third
        ),
        remaining_quintic_proxy_entries=remaining_quintic,
        remaining_quintic_local_proxy=local_proxy,
        proxy_reserve_after_declared_allowance=(
            local_proxy.margin_to_one_third - RESERVE_TARGET
        ),
    )


def artifact_text(result: FixedPairAdjacentRowContraction) -> str:
    payload = {
        "schema": "round4_q64_fixed_pair_adjacent_row_contraction_v1",
        "result": asdict(result),
        "evidence_label": (
            "arbitrary-diagonal complete-row Schur-feature theorem using "
            "fixed-pair cubic completion counts and inherited record-specific "
            "M_35 pair-tail bounds; floating q64 Perron insertion; one batch only"
        ),
        "remaining_open": (
            "624 balanced entries plus interval certification and the adaptive lift"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic fixed-pair adjacent contraction",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_fixed_pair_adjacent_row_contraction.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 fixed-pair adjacent row contraction: "
        f"entries={result.closed_entries},"
        f"coefficient={result.coefficient:.12g},"
        f"proved_entries={result.total_proved_entries},"
        f"remaining={result.remaining_open_entries},"
        f"total={result.fixed_pair_inserted.total:.12g},"
        f"margin={result.fixed_pair_inserted.margin_to_one_third:.12g},"
        f"proxy_total={result.remaining_quintic_local_proxy.total:.12g},"
        "proxy_reserve_after_allowance="
        f"{result.proxy_reserve_after_declared_allowance:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
