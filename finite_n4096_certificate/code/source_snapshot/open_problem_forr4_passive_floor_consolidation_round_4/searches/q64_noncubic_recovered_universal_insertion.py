#!/usr/bin/env python3
"""Close the q64 noncubic class and recover one universal residual class."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent
    / "open_problem_forr4_passive_floor_consolidation_round_3"
    / "searches"
)
sys.path.insert(0, str(ROUND3_SEARCHES))

from high_degree_record_incidence_frontier import (  # noqa: E402
    sector_coefficient_bound,
)
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
from q64_remaining_class_gates import partition_remaining  # noqa: E402
from q64_shared_quintic_row_chain_insertion import (  # noqa: E402
    inserted_coefficients as shared_quintic_inserted_coefficients,
)
from q64_universal_septimic_insertion import (  # noqa: E402
    UNIVERSAL_GRAM_COEFFICIENT,
)


NONCUBIC_CLASS = "noncubic_profile"
RECOVERED_CLASS = "two_split_cubics_one_split_higher"
REMAINING_CLASS_LABELS = (
    "higher_split_only_in_cubic_profile",
    "one_split_cubic_no_split_higher",
)
REMAINING_CLASS_RESERVE_GATES = (
    0.14034303056468567,
    0.3448872174126536,
)


@dataclass(frozen=True)
class NoncubicRecoveredUniversalInsertion:
    order: int
    dimension: int
    sign_modes: int
    injective_character_product_factor: float
    degree_seven_entries: int
    degree_seven_orbits: int
    degree_seven_minimum_coefficient: float
    degree_seven_maximum_coefficient: float
    universal_noncubic_entries: int
    noncubic_entries: int
    recovered_universal_entries: int
    recovered_universal_coefficient: float
    previous_proved_entries: int
    post_noncubic_proved_entries: int
    total_proved_entries: int
    remaining_open_entries: int
    previous_routing: OptimizedLedger
    post_noncubic_routing: OptimizedLedger
    final_routing: OptimizedLedger
    noncubic_routing_improvement: float
    recovered_universal_routing_cost: float
    net_routing_change: float
    reserve_after_declared_allowance: float
    adaptive_multiplier_cap_retaining_allowance: float
    remaining_class_labels: tuple[str, ...]
    remaining_class_counts: tuple[int, ...]
    remaining_class_frozen_targets: tuple[float, ...]
    remaining_class_reserve_gates: tuple[float, ...]


def character_product_factor(order: int = ORDER) -> float:
    """Worst product of the two degree-seven endpoint residuals."""

    q = order
    if q < 4 or q & (q - 1):
        raise ValueError(("power-of-two order at least four required", q))
    return 3 / ((q - 1) * (q - 3))


def noncubic_entries() -> tuple[ProfileSplit, ...]:
    return partition_remaining()[NONCUBIC_CLASS]


def degree_seven_entries() -> tuple[ProfileSplit, ...]:
    return tuple(entry for entry in noncubic_entries() if 7 in entry[0])


def universal_noncubic_entries() -> tuple[ProfileSplit, ...]:
    return tuple(entry for entry in noncubic_entries() if 7 not in entry[0])


def recovered_universal_entries() -> tuple[ProfileSplit, ...]:
    return partition_remaining()[RECOVERED_CLASS]


def degree_seven_coefficient(
    entry: ProfileSplit,
    order: int = ORDER,
) -> float:
    """Rank/incidence coefficient with the joint endpoint improvement."""

    profile, _ = entry
    if 7 not in profile or profile.count(1) != 3:
        raise ValueError(("not a single middle septimic entry", entry))
    # Every link record is forced to one. The inherited sector estimate uses
    # the safe entry cap q^-3. The endpoint-product lemma improves that cap
    # by exactly ``character_product_factor`` without changing its geometry.
    return (
        sector_coefficient_bound(order, *entry, (1, 1, 1))
        * character_product_factor(order)
    )


def noncubic_coefficient_map() -> dict[ProfileSplit, float]:
    result = {
        entry: degree_seven_coefficient(entry)
        for entry in degree_seven_entries()
    }
    result.update(
        {
            entry: UNIVERSAL_GRAM_COEFFICIENT
            for entry in universal_noncubic_entries()
        }
    )
    return result


def post_noncubic_coefficients() -> dict[ProfileSplit, float]:
    result = shared_quintic_inserted_coefficients()
    result.update(noncubic_coefficient_map())
    return result


def inserted_coefficients() -> dict[ProfileSplit, float]:
    result = post_noncubic_coefficients()
    for entry in recovered_universal_entries():
        result[entry] = UNIVERSAL_GRAM_COEFFICIENT
    return result


def diagnostic() -> NoncubicRecoveredUniversalInsertion:
    degree_seven = degree_seven_entries()
    degree_seven_coefficients = tuple(
        degree_seven_coefficient(entry) for entry in degree_seven
    )
    universal_noncubic = universal_noncubic_entries()
    noncubic = noncubic_entries()
    recovered = recovered_universal_entries()
    previous = optimize(
        mapped_coefficients=shared_quintic_inserted_coefficients()
    )
    post_noncubic = optimize(mapped_coefficients=post_noncubic_coefficients())
    final = optimize(mapped_coefficients=inserted_coefficients())
    previous_proved = 428
    remaining_partition = partition_remaining()
    remaining_classes = tuple(
        remaining_partition[label] for label in REMAINING_CLASS_LABELS
    )
    final_coefficients = inserted_coefficients()
    return NoncubicRecoveredUniversalInsertion(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        injective_character_product_factor=character_product_factor(),
        degree_seven_entries=len(degree_seven),
        degree_seven_orbits=len(degree_seven) // 4,
        degree_seven_minimum_coefficient=min(degree_seven_coefficients),
        degree_seven_maximum_coefficient=max(degree_seven_coefficients),
        universal_noncubic_entries=len(universal_noncubic),
        noncubic_entries=len(noncubic),
        recovered_universal_entries=len(recovered),
        recovered_universal_coefficient=UNIVERSAL_GRAM_COEFFICIENT,
        previous_proved_entries=previous_proved,
        post_noncubic_proved_entries=previous_proved + len(noncubic),
        total_proved_entries=previous_proved + len(noncubic) + len(recovered),
        remaining_open_entries=(
            888 - previous_proved - len(noncubic) - len(recovered)
        ),
        previous_routing=previous,
        post_noncubic_routing=post_noncubic,
        final_routing=final,
        noncubic_routing_improvement=previous.total - post_noncubic.total,
        recovered_universal_routing_cost=final.total - post_noncubic.total,
        net_routing_change=final.total - previous.total,
        reserve_after_declared_allowance=(
            final.margin_to_one_third - RESERVE_TARGET
        ),
        adaptive_multiplier_cap_retaining_allowance=(
            (THRESHOLD - RESERVE_TARGET) / final.total
        ),
        remaining_class_labels=REMAINING_CLASS_LABELS,
        remaining_class_counts=tuple(
            len(entries) for entries in remaining_classes
        ),
        remaining_class_frozen_targets=tuple(
            final_coefficients[entries[0]] for entries in remaining_classes
        ),
        remaining_class_reserve_gates=REMAINING_CLASS_RESERVE_GATES,
    )


def artifact_text(result: NoncubicRecoveredUniversalInsertion) -> str:
    payload = {
        "schema": "round4_q64_noncubic_recovered_universal_insertion_v1",
        "result": asdict(result),
        "evidence_label": (
            "mixed result: the specialized arbitrary-diagonal theorem for 16 "
            "middle-septimic entries survives; 124 noncubic and 96 recovered "
            "universal entries are quarantined because the cross-Gram step "
            "omits distinctness masks"
        ),
        "adaptive_requirement": (
            "the displayed reserve and multiplier are conditional on the "
            "frozen targets for 224 remaining entries; they are not an "
            "adaptive recurrence or passive theorem"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic noncubic/recovered insertion artifact",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = (
            ROOT
            / "artifacts"
            / "q64_noncubic_recovered_universal_insertion.json"
        )
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 noncubic/recovered universal insertion: "
        f"noncubic={result.noncubic_entries},"
        f"degree7_range={result.degree_seven_minimum_coefficient:.12g}/"
        f"{result.degree_seven_maximum_coefficient:.12g},"
        f"recovered={result.recovered_universal_entries},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.final_routing.total:.12g},"
        f"margin={result.final_routing.margin_to_one_third:.12g},"
        f"reserve={result.reserve_after_declared_allowance:.12g},"
        f"remaining_open={result.remaining_open_entries},"
        "status=mixed_16_proved_220_quarantined_cumulative_values_withdrawn"
    )


if __name__ == "__main__":
    main()
