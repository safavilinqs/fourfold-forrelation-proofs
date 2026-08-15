#!/usr/bin/env python3
"""Exact-order routing gate for a possible N=4096 paper result.

The inherited occupation certificate was written with module-level q=32
constants even though its coefficient formulas are valid for general
power-of-two order.  This script evaluates those same formulas at q=64,
retains all already proved degree-four through known-high-sector bounds, and
asks how much common coefficient the still-open degree-ten/twelve profiles
could tolerate.

It also evaluates one explicit two-tier *proof target*: open profiles that
contain a cubic block are assigned the q=64 cubic fixed-pair slice scale,
while the remaining open profiles are assigned the deliberately loose value
1/2.  Neither assignment is asserted as a theorem here.  The output is a
floating routing certificate for choosing the next contraction and size.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from json import dumps
from math import sqrt
from pathlib import Path
import sys
from typing import Iterator

from scipy.optimize import brentq, minimize_scalar


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent / "open_problem_forr4_passive_floor_consolidation_round_3" / "searches"
)
sys.path.insert(0, str(ROUND3_SEARCHES))

import degree_six_joint_occupation_optimization as degree_six  # noqa: E402
import double_endpoint_occupation_optimization as double_endpoint  # noqa: E402
import occupation_compatible_sector_optimization as occupation  # noqa: E402
from attenuation_promise_concentration import (  # noqa: E402
    extended_euclidean_promise_concentration,
)


ORDER = 64
DIMENSION = ORDER * ORDER
MODES = 4 * DIMENSION
THRESHOLD = 1 / 3
RESERVE_TARGET = 1e-3
BETA_BOUNDS = (2 ** (-0.5) + 1e-5, 0.82)

Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class OptimizedLedger:
    total: float
    beta: float
    perron_upper: float
    promise_loss: float
    margin_to_one_third: float


@dataclass(frozen=True)
class Q64TargetGate:
    order: int
    dimension: int
    sign_modes: int
    balanced_open_entries: int
    cubic_profile_entries: int
    noncubic_profile_entries: int
    block_coherent_entries: int
    internally_split_cubic_entries: int
    cubic_profile_high_internal_entries: int
    noncubic_high_internal_entries: int
    common_one_over_q: OptimizedLedger
    common_one_over_sqrt_q: OptimizedLedger
    common_threshold_coefficient: float
    common_threshold_times_q: float
    common_reserve_coefficient: float
    common_reserve_times_q: float
    cubic_slice_target: float
    noncubic_target: float
    two_tier_target: OptimizedLedger


@contextmanager
def order_context(order: int) -> Iterator[None]:
    """Temporarily evaluate inherited formulas at a different order."""

    modules = (occupation, degree_six, double_endpoint)
    previous = tuple(module.ORDER for module in modules)
    try:
        for module in modules:
            module.ORDER = order
        yield
    finally:
        for module, value in zip(modules, previous, strict=True):
            module.ORDER = value


def open_profiles() -> tuple[Profile, ...]:
    return tuple(
        profile
        for profile in occupation.HIGH_DEGREE_PROFILES
        if profile not in occupation.KNOWN_HIGH_DEGREE_PROFILES
    )


def balanced_open_entries() -> tuple[ProfileSplit, ...]:
    return tuple(
        (profile, split)
        for profile in open_profiles()
        for split in occupation.profile_splits(profile)
        if 2 * sum(split) == sum(profile)
    )


def optimize(
    common_coefficient: float | None = None,
    mapped_coefficients: dict[ProfileSplit, float] | None = None,
) -> OptimizedLedger:
    if (common_coefficient is None) == (mapped_coefficients is None):
        raise ValueError("specify exactly one high-sector coefficient model")

    with order_context(ORDER):

        def total(beta: float) -> float:
            ledger = occupation.certificate(
                beta=beta,
                high_degree_coefficient=common_coefficient,
                profile_split_coefficients=mapped_coefficients,
            )
            promise = extended_euclidean_promise_concentration(DIMENSION, beta)
            return ledger.supporting_upper + promise.two_hypothesis_loss

        optimum = minimize_scalar(
            total,
            bounds=BETA_BOUNDS,
            method="bounded",
            options={"xatol": 2e-10},
        )
        beta = float(optimum.x)
        ledger = occupation.certificate(
            beta=beta,
            high_degree_coefficient=common_coefficient,
            profile_split_coefficients=mapped_coefficients,
        )
        promise = extended_euclidean_promise_concentration(DIMENSION, beta)
        result = ledger.supporting_upper + promise.two_hypothesis_loss
    return OptimizedLedger(
        total=result,
        beta=beta,
        perron_upper=ledger.supporting_upper,
        promise_loss=promise.two_hypothesis_loss,
        margin_to_one_third=THRESHOLD - result,
    )


def common_gate(target_total: float) -> float:
    cache: dict[float, float] = {}

    def residual(value: float) -> float:
        if value not in cache:
            cache[value] = optimize(common_coefficient=value).total
        return cache[value] - target_total

    upper = 0.25
    while residual(upper) < 0:
        upper *= 2
    return float(brentq(residual, 0.0, upper, xtol=2e-9))


def two_tier_coefficients(
    cubic_value: float, noncubic_value: float
) -> dict[ProfileSplit, float]:
    return {
        (profile, split): cubic_value if 3 in profile else noncubic_value
        for profile in open_profiles()
        for split in occupation.profile_splits(profile)
    }


def contraction_class(entry: ProfileSplit) -> str:
    profile, split = entry
    if all(
        selected in (0, degree) for degree, selected in zip(profile, split, strict=True)
    ):
        return "block_coherent"
    if any(
        degree == 3 and selected in (1, 2)
        for degree, selected in zip(profile, split, strict=True)
    ):
        return "internally_split_cubic"
    if 3 in profile:
        return "cubic_profile_high_internal"
    return "noncubic_high_internal"


def diagnostic() -> Q64TargetGate:
    entries = balanced_open_entries()
    cubic_entries = tuple(entry for entry in entries if 3 in entry[0])
    noncubic_entries = tuple(entry for entry in entries if 3 not in entry[0])
    classes = {
        label: 0
        for label in (
            "block_coherent",
            "internally_split_cubic",
            "cubic_profile_high_internal",
            "noncubic_high_internal",
        )
    }
    for entry in entries:
        classes[contraction_class(entry)] += 1
    with order_context(ORDER):
        cubic_slice = sqrt(occupation.endpoint_singleton_slice_energies(ORDER)[2])
    noncubic_target = 0.5
    two_tier = two_tier_coefficients(cubic_slice, noncubic_target)
    threshold_gate = common_gate(THRESHOLD)
    reserve_gate = common_gate(THRESHOLD - RESERVE_TARGET)
    return Q64TargetGate(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        balanced_open_entries=len(entries),
        cubic_profile_entries=len(cubic_entries),
        noncubic_profile_entries=len(noncubic_entries),
        block_coherent_entries=classes["block_coherent"],
        internally_split_cubic_entries=classes["internally_split_cubic"],
        cubic_profile_high_internal_entries=classes["cubic_profile_high_internal"],
        noncubic_high_internal_entries=classes["noncubic_high_internal"],
        common_one_over_q=optimize(common_coefficient=1 / ORDER),
        common_one_over_sqrt_q=optimize(common_coefficient=1 / sqrt(ORDER)),
        common_threshold_coefficient=threshold_gate,
        common_threshold_times_q=threshold_gate * ORDER,
        common_reserve_coefficient=reserve_gate,
        common_reserve_times_q=reserve_gate * ORDER,
        cubic_slice_target=cubic_slice,
        noncubic_target=noncubic_target,
        two_tier_target=optimize(mapped_coefficients=two_tier),
    )


def artifact_text(result: Q64TargetGate) -> str:
    payload = {
        "schema": "round4_q64_paper_target_gate_v1",
        "result": asdict(result),
        "evidence_label": (
            "inherited proved lower-sector formulas evaluated at q=64; "
            "floating Perron and attenuation optimization; high-sector "
            "coefficients are routing targets, not proved bounds; one batch only"
        ),
        "retain_condition": (
            "prove the displayed common gate or the sufficient two-tier "
            "envelope, then intervalize and lift to adaptive protocols"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic q64 routing result under artifacts/",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_paper_target_gate.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 paper target gate: "
        f"N={result.dimension},"
        f"M={result.sign_modes},"
        f"open_entries={result.balanced_open_entries},"
        f"common_1_over_q_total={result.common_one_over_q.total:.12g},"
        f"common_1_over_sqrt_q_total={result.common_one_over_sqrt_q.total:.12g},"
        f"common_gate={result.common_threshold_coefficient:.12g},"
        f"reserve_gate={result.common_reserve_coefficient:.12g},"
        f"cubic_target={result.cubic_slice_target:.12g},"
        f"noncubic_target={result.noncubic_target:.12g},"
        f"two_tier_total={result.two_tier_target.total:.12g},"
        f"two_tier_margin={result.two_tier_target.margin_to_one_third:.12g},"
        "status=routing_target_not_theorem"
    )


if __name__ == "__main__":
    main()
