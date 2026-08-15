#!/usr/bin/env python3
"""Finite-size scalar-ledger screen for the quadratic-bent replacement.

The quadratic-bent plant has exact promise separation, so its natural first
screen uses ``beta=1`` and no conditioning penalty. This program asks how
small a *complete* family of higher-sector coefficients would have to be at
``N=1024`` before that promise advantage yields passive dose greater than
six in the current occupation-compatible scalar architecture.

This is deliberately an optimistic screen. Unknown higher-sector
coefficients are set to zero when computing the floor and to a common value
when computing the gate. It does not assert that the quadratic-bent plant
has those coefficients, and it is not a passive theorem.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from itertools import product
from json import dumps
from pathlib import Path
import sys

from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent / "open_problem_forr4_passive_floor_consolidation_round_3" / "searches"
)
sys.path.insert(0, str(ROUND3_SEARCHES))

from double_endpoint_occupation_optimization import (  # noqa: E402
    occupation_states,
)
from occupation_compatible_sector_optimization import (  # noqa: E402
    certificate,
    multiplicity,
    paired_state,
    profile_splits,
)


DIMENSION = 1024
THRESHOLD = 1 / 3
HIGHER_DEGREES = (6, 8, 10, 12)
ODD_BLOCK_DEGREES = (1, 3, 5, 7, 9, 11)

Profile = tuple[int, ...]
Split = tuple[int, ...]
Entry = tuple[Profile, Split]


@dataclass(frozen=True)
class ReplacementScreen:
    dimension: int
    beta: float
    promise_loss: float
    nonminimal_profiles: int
    raw_profile_splits: int
    compatible_profile_splits: int
    compatible_by_degree: dict[int, int]
    optimistic_zero_higher_floor: float
    optimistic_floor_margin: float
    common_coefficient_gate: float
    common_gate_times_dimension: float
    half_over_dimension_total: float
    half_over_dimension_margin: float
    one_over_dimension_total: float
    one_over_dimension_overshoot: float
    degree_only_gates: dict[int, float]
    known_endpoint_five_one_coefficient: float
    endpoint_to_common_gate_ratio: float
    decision: str
    evidence_label: str


def profiles() -> tuple[Profile, ...]:
    return tuple(
        profile
        for profile in product(ODD_BLOCK_DEGREES, repeat=4)
        if sum(profile) in HIGHER_DEGREES
    )


def raw_entries() -> tuple[Entry, ...]:
    return tuple(
        (profile, split) for profile in profiles() for split in profile_splits(profile)
    )


def compatible_entries() -> tuple[Entry, ...]:
    states = occupation_states()
    state_set = set(states)
    result: list[Entry] = []
    for profile, split in raw_entries():
        complement = tuple(
            degree - selected for degree, selected in zip(profile, split, strict=True)
        )
        compatible = False
        for state in states:
            if any(
                occupation < selected
                for occupation, selected in zip(state, split, strict=True)
            ):
                continue
            partner = paired_state(state, profile, split)
            if partner not in state_set:
                continue
            if multiplicity(state, split) and multiplicity(partner, complement):
                compatible = True
                break
        if compatible:
            result.append((profile, split))
    return tuple(result)


def scalar_total(entries: tuple[Entry, ...], coefficients: dict[int, float]) -> float:
    mapped = {entry: coefficients.get(sum(entry[0]), 0.0) for entry in entries}
    return certificate(
        beta=1.0,
        profile_split_coefficients=mapped,
    ).supporting_upper


def screen() -> ReplacementScreen:
    profile_list = profiles()
    entry_list = raw_entries()
    compatible = compatible_entries()
    by_degree = {
        degree: sum(sum(profile) == degree for profile, _ in compatible)
        for degree in HIGHER_DEGREES
    }

    floor = scalar_total(compatible, {})

    def common_total(coefficient: float) -> float:
        return scalar_total(
            compatible,
            {degree: coefficient for degree in HIGHER_DEGREES},
        )

    common_gate = brentq(
        lambda coefficient: common_total(coefficient) - THRESHOLD,
        0.0,
        1 / 100,
        xtol=1e-14,
    )
    degree_gates = {
        degree: brentq(
            lambda coefficient, active=degree: (
                scalar_total(compatible, {active: coefficient}) - THRESHOLD
            ),
            0.0,
            1 / 10,
            xtol=1e-14,
        )
        for degree in HIGHER_DEGREES
    }
    half_total = common_total(1 / (2 * DIMENSION))
    one_total = common_total(1 / DIMENSION)
    endpoint = 2 / (DIMENSION - 2)
    return ReplacementScreen(
        dimension=DIMENSION,
        beta=1.0,
        promise_loss=0.0,
        nonminimal_profiles=len(profile_list),
        raw_profile_splits=len(entry_list),
        compatible_profile_splits=len(compatible),
        compatible_by_degree=by_degree,
        optimistic_zero_higher_floor=floor,
        optimistic_floor_margin=THRESHOLD - floor,
        common_coefficient_gate=common_gate,
        common_gate_times_dimension=common_gate * DIMENSION,
        half_over_dimension_total=half_total,
        half_over_dimension_margin=THRESHOLD - half_total,
        one_over_dimension_total=one_total,
        one_over_dimension_overshoot=one_total - THRESHOLD,
        degree_only_gates=degree_gates,
        known_endpoint_five_one_coefficient=endpoint,
        endpoint_to_common_gate_ratio=endpoint / common_gate,
        decision="not_promoted_current_signed_permutation_remains_lead",
        evidence_label=(
            "exact scalar-ledger sensitivity calculation with hypothetical "
            "common higher-sector coefficients; replacement screen, not "
            "quadratic-bent coefficient theorem"
        ),
    )


def artifact_text(result: ReplacementScreen) -> str:
    return (
        dumps(
            {
                "schema": "round4_quadratic_bent_replacement_screen_v1",
                **asdict(result),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic replacement-screen JSON artifact",
    )
    arguments = parser.parse_args()
    result = screen()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "quadratic_bent_replacement_screen.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "quadratic-bent replacement screen: "
        f"N={result.dimension},"
        f"profiles={result.nonminimal_profiles},"
        f"raw_entries={result.raw_profile_splits},"
        f"compatible_entries={result.compatible_profile_splits},"
        f"floor={result.optimistic_zero_higher_floor:.12g},"
        f"floor_margin={result.optimistic_floor_margin:.12g},"
        f"common_gate={result.common_coefficient_gate:.12g},"
        f"gate_times_N={result.common_gate_times_dimension:.12g},"
        f"half_over_N_total={result.half_over_dimension_total:.12g},"
        f"half_over_N_margin={result.half_over_dimension_margin:.12g},"
        f"one_over_N_total={result.one_over_dimension_total:.12g},"
        f"one_over_N_overshoot={result.one_over_dimension_overshoot:.12g},"
        f"endpoint_5_1={result.known_endpoint_five_one_coefficient:.12g},"
        f"decision={result.decision}"
    )


if __name__ == "__main__":
    main()
