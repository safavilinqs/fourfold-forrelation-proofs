#!/usr/bin/env python3
"""Record the repair audit triggered by the adjacent direct contraction.

Both certified physical witness orbits remain fixed, and the leading adjacent
split orbit is charged at its exact combined local coefficient.  For each
already accepted profile family, uniformly scale its proved fixed-split
coefficients, reoptimize attenuation using the two-split promise theorem, and
solve for the smallest reduction that restores the one-third threshold.

Historically this calculation ranked accepted-sector repairs.  The selected
adjacent double-cubic repair is now implemented in the live coefficient
ledger, so the current main routine first verifies that no further repair is
needed and reports the realized margin.
"""

from __future__ import annotations

from dataclasses import dataclass

from scipy.optimize import brentq, minimize_scalar

from attenuation_promise_concentration import (
    extended_euclidean_promise_concentration,
)
from occupation_compatible_sector_optimization import (
    SEPARATED_CUBIC_QUINTIC_PROFILES,
    SEPARATED_QUINTIC_CUBIC_PROFILES,
    TRIPLE_CUBIC_PROFILES,
    certificate,
    coefficient,
    profile_splits,
)
from repaired_open_profile_budget import (
    ADJACENT_SPLIT_ENTRIES,
    DIMENSION,
    THRESHOLD,
    hybrid_forced_coefficients,
    split_orbits,
)
from adjacent_cubic_quintic_orbit_witness import (
    horizontal_adjacent_slice_certificate,
)


@dataclass(frozen=True)
class AdjacentRepairTarget:
    name: str
    profiles: tuple[tuple[int, ...], ...]
    required_reduction: float
    optimal_beta: float
    threshold_scale: float


@dataclass(frozen=True)
class AdjacentSplitRepairTarget:
    entries: tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]
    removal_gain: float
    linear_required_reduction: float


GROUPS = {
    "adjacent_double_cubic": ((1, 1, 3, 3), (3, 3, 1, 1)),
    "endpoint_cubic_degree_six": ((3, 1, 1, 1), (1, 1, 1, 3)),
    "endpoint_quintic_degree_eight": ((5, 1, 1, 1), (1, 1, 1, 5)),
    "all_triple_cubic": TRIPLE_CUBIC_PROFILES,
    "central_double_cubic": ((1, 3, 3, 1),),
    "double_endpoint": ((3, 1, 1, 3),),
    "separated_cubic_quintic": (
        SEPARATED_QUINTIC_CUBIC_PROFILES
        + SEPARATED_CUBIC_QUINTIC_PROFILES
    ),
}


def adjacent_forced_coefficients() -> dict[
    tuple[tuple[int, ...], tuple[int, ...]], float
]:
    result = hybrid_forced_coefficients()
    combined = horizontal_adjacent_slice_certificate(32).combined_coefficient
    for entry in ADJACENT_SPLIT_ENTRIES:
        result.setdefault(entry, combined)
    return result


def optimized_total(
    coefficients: dict[tuple[tuple[int, ...], tuple[int, ...]], float],
) -> tuple[float, float]:
    def total(beta: float) -> float:
        ledger = certificate(
            beta=beta,
            profile_split_coefficients=coefficients,
        )
        promise = extended_euclidean_promise_concentration(DIMENSION, beta)
        return ledger.supporting_upper + promise.two_hypothesis_loss

    result = minimize_scalar(
        total,
        bounds=(0.75, 0.81),
        method="bounded",
        options={"xatol": 1e-12},
    )
    return float(result.fun), float(result.x)


def fixed_total(
    beta: float,
    coefficients: dict[tuple[tuple[int, ...], tuple[int, ...]], float],
) -> float:
    ledger = certificate(
        beta=beta,
        profile_split_coefficients=coefficients,
    )
    promise = extended_euclidean_promise_concentration(DIMENSION, beta)
    return ledger.supporting_upper + promise.two_hypothesis_loss


def repair_target(
    name: str,
    profiles: tuple[tuple[int, ...], ...],
    forced: dict[tuple[tuple[int, ...], tuple[int, ...]], float],
) -> AdjacentRepairTarget:
    proved = {
        (profile, split): coefficient(profile, split, True, None)
        for profile in profiles
        for split in profile_splits(profile)
    }

    def mapped(scale: float):
        result = dict(forced)
        result.update({key: scale * value for key, value in proved.items()})
        return result

    if optimized_total(mapped(0))[0] >= THRESHOLD:
        raise AssertionError(("profile family cannot repair budget", name))
    scale = brentq(
        lambda value: optimized_total(mapped(value))[0] - THRESHOLD,
        0,
        1,
        xtol=2e-9,
    )
    _, beta = optimized_total(mapped(scale))
    return AdjacentRepairTarget(
        name=name,
        profiles=profiles,
        required_reduction=1 - scale,
        optimal_beta=beta,
        threshold_scale=scale,
    )


def adjacent_double_cubic_split_targets(
    forced: dict[tuple[tuple[int, ...], tuple[int, ...]], float],
    beta: float,
) -> tuple[AdjacentSplitRepairTarget, ...]:
    profiles = GROUPS["adjacent_double_cubic"]
    proved = {
        (profile, split): coefficient(profile, split, True, None)
        for profile in profiles
        for split in profile_splits(profile)
    }
    baseline_coefficients = dict(forced)
    baseline_coefficients.update(proved)
    baseline = fixed_total(beta, baseline_coefficients)
    overshoot = baseline - THRESHOLD
    result = []
    for entries in split_orbits(profiles[0]):
        modified = dict(baseline_coefficients)
        for entry in entries:
            modified[entry] = 0.0
        gain = baseline - fixed_total(beta, modified)
        if gain <= 1e-14:
            continue
        result.append(
            AdjacentSplitRepairTarget(
                entries=entries,
                removal_gain=gain,
                linear_required_reduction=overshoot / gain,
            )
        )
    return tuple(
        sorted(result, key=lambda item: item.removal_gain, reverse=True)
    )


def main() -> None:
    forced = adjacent_forced_coefficients()
    baseline, beta = optimized_total(forced)
    print(
        "adjacent combined diagnostic after accepted-sector repair: "
        f"beta={beta:.12g},total={baseline:.12g},"
        f"slack={THRESHOLD - baseline:.12g}"
    )
    if baseline < THRESHOLD:
        return
    targets = sorted(
        (
            repair_target(name, profiles, forced)
            for name, profiles in GROUPS.items()
        ),
        key=lambda result: result.required_reduction,
    )
    for target in targets:
        print(
            "adjacent accepted-sector repair target: "
            f"name={target.name},"
            f"required_reduction={target.required_reduction:.12g},"
            f"threshold_scale={target.threshold_scale:.12g},"
            f"optimal_beta={target.optimal_beta:.12g}"
        )
    for target in adjacent_double_cubic_split_targets(forced, beta)[:8]:
        print(
            "adjacent double-cubic split repair target: "
            f"entries={target.entries},"
            f"removal_gain={target.removal_gain:.12g},"
            f"linear_required_reduction="
            f"{target.linear_required_reduction:.12g}"
        )


if __name__ == "__main__":
    main()
