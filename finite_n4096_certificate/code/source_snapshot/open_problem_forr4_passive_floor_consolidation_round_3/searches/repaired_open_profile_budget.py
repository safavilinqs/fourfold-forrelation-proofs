#!/usr/bin/env python3
"""Route-selection ledgers for the open high-degree profiles.

The explicit opposite-endpoint mixed-orbit cuts are held at their proved
q=32 witness value.  Every other open split is initially zero.  This script
then measures the historical promise-repair thresholds and the superseding
coarse completion target after the chained accepted-sector repair:

1. the common coefficient that all remaining splits could tolerate;
2. the first-order cost of each profile reversal orbit; and
3. the first-order cost of each split orbit in the leading profile.
4. the all-open ledger obtained by assigning 1/q to every unforced split.

These are route-selection diagnostics.  A threshold is not a proved tester
coefficient.
"""

from __future__ import annotations

from dataclasses import dataclass

from scipy.optimize import brentq, minimize_scalar

from attenuation_promise_concentration import (
    euclidean_promise_concentration,
    extended_euclidean_promise_concentration,
    hybrid_euclidean_promise_concentration,
)
from adjacent_cubic_quintic_orbit_witness import (
    horizontal_adjacent_slice_certificate,
)
from occupation_compatible_sector_optimization import (
    HIGH_DEGREE_PROFILES,
    KNOWN_HIGH_DEGREE_PROFILES,
    certificate,
    profile_splits,
)
from opposite_endpoint_vertical_mixture_witness import (
    forced_split_coefficients,
)


DIMENSION = 1024
WITNESS_COEFFICIENT = 0.039593955294628
REPAIRED_BETA = 0.779512326135891
CERTIFIED_LEADING_COEFFICIENT = 0.014281024204693649
HYBRID_REPAIRED_BETA = 0.7797679764624451
EXTENDED_REPAIRED_BETA = 0.7798687433309669
THRESHOLD = 1 / 3


Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]
LEADING_SPLIT_ENTRIES: tuple[ProfileSplit, ...] = (
    ((3, 1, 1, 5), (1, 0, 1, 3)),
    ((3, 1, 1, 5), (2, 1, 0, 2)),
    ((5, 1, 1, 3), (2, 0, 1, 2)),
    ((5, 1, 1, 3), (3, 1, 0, 1)),
)
ADJACENT_SPLIT_ENTRIES: tuple[ProfileSplit, ...] = (
    ((1, 1, 3, 5), (0, 1, 2, 2)),
    ((1, 1, 3, 5), (1, 0, 1, 3)),
    ((5, 3, 1, 1), (2, 2, 1, 0)),
    ((5, 3, 1, 1), (3, 1, 0, 1)),
)


@dataclass(frozen=True)
class RepairedOpenBudget:
    beta: float
    baseline_total: float
    baseline_slack: float
    common_coefficient_threshold: float


@dataclass(frozen=True)
class OrbitSensitivity:
    profiles: tuple[Profile, ...]
    slope: float
    linear_threshold: float


@dataclass(frozen=True)
class SplitOrbitSensitivity:
    entries: tuple[ProfileSplit, ...]
    slope: float
    linear_threshold: float


@dataclass(frozen=True)
class LeadingScalarObstruction:
    coefficient: float
    optimal_beta: float
    optimized_total: float
    threshold_overshoot: float


@dataclass(frozen=True)
class LeadingHybridRepair:
    coefficient: float
    optimal_beta: float
    optimized_total: float
    threshold_slack: float


@dataclass(frozen=True)
class ExtendedAdjacentDiagnostic:
    scalar_threshold: float
    threshold_beta: float
    record_one_coefficient: float
    record_three_coefficient: float
    combined_slice_coefficient: float
    record_one_total: float
    record_one_slack: float
    combined_total: float
    combined_overshoot: float


@dataclass(frozen=True)
class CoarseOpenCompletionTarget:
    """Diagnostic obtained by charging every unforced open split at 1/q.

    This is a proof target, not a theorem that the open coefficients obey
    the stated bound.
    """

    coefficient: float
    optimal_beta: float
    optimized_total: float
    threshold_slack: float
    forced_entries: int
    coarse_entries: int


def open_profiles() -> tuple[Profile, ...]:
    return tuple(
        profile
        for profile in HIGH_DEGREE_PROFILES
        if profile not in KNOWN_HIGH_DEGREE_PROFILES
    )


def reversal_orbits(profiles: tuple[Profile, ...]) -> tuple[tuple[Profile, ...], ...]:
    available = set(profiles)
    seen: set[Profile] = set()
    result = []
    for profile in profiles:
        if profile in seen:
            continue
        orbit = tuple(
            sorted(
                {
                    profile,
                    tuple(reversed(profile)),
                }
                & available
            )
        )
        seen.update(orbit)
        result.append(orbit)
    return tuple(result)


def forced_coefficients() -> dict[ProfileSplit, float]:
    return forced_split_coefficients(WITNESS_COEFFICIENT)


def mapped_coefficients(
    profiles: tuple[Profile, ...],
    value: float,
) -> dict[ProfileSplit, float]:
    result = forced_coefficients()
    for profile in profiles:
        for split in profile_splits(profile):
            result.setdefault((profile, split), value)
    return result


def fixed_beta_total(coefficients: dict[ProfileSplit, float]) -> float:
    ledger = certificate(
        beta=REPAIRED_BETA,
        profile_split_coefficients=coefficients,
    )
    promise = euclidean_promise_concentration(DIMENSION, REPAIRED_BETA)
    return ledger.supporting_upper + promise.two_hypothesis_loss


def baseline_total() -> float:
    return fixed_beta_total(forced_coefficients())


def common_threshold(profiles: tuple[Profile, ...]) -> float:
    def residual(value: float) -> float:
        return fixed_beta_total(mapped_coefficients(profiles, value)) - THRESHOLD

    upper = 1 / 10
    if residual(0) >= 0 or residual(upper) <= 0:
        raise AssertionError(("common-threshold bracket", residual(0), residual(upper)))
    return float(brentq(residual, 0, upper, xtol=1e-13))


def repaired_open_budget() -> RepairedOpenBudget:
    profiles = open_profiles()
    baseline = baseline_total()
    return RepairedOpenBudget(
        beta=REPAIRED_BETA,
        baseline_total=baseline,
        baseline_slack=THRESHOLD - baseline,
        common_coefficient_threshold=common_threshold(profiles),
    )


def orbit_sensitivities(step: float = 1e-5) -> tuple[OrbitSensitivity, ...]:
    baseline = baseline_total()
    result = []
    for orbit in reversal_orbits(open_profiles()):
        slope = (
            fixed_beta_total(mapped_coefficients(orbit, step)) - baseline
        ) / step
        result.append(
            OrbitSensitivity(
                profiles=orbit,
                slope=slope,
                linear_threshold=(
                    (THRESHOLD - baseline) / slope if slope > 1e-12 else float("inf")
                ),
            )
        )
    return tuple(sorted(result, key=lambda item: item.slope, reverse=True))


def split_orbits(profile: Profile) -> tuple[tuple[ProfileSplit, ...], ...]:
    reverse = tuple(reversed(profile))
    seen: set[ProfileSplit] = set()
    result = []
    for split in profile_splits(profile):
        complement = tuple(
            degree - selected
            for degree, selected in zip(profile, split, strict=True)
        )
        orbit = tuple(
            sorted(
                {
                    (profile, split),
                    (profile, complement),
                    (reverse, tuple(reversed(split))),
                    (reverse, tuple(reversed(complement))),
                }
            )
        )
        if any(entry in seen for entry in orbit):
            continue
        seen.update(orbit)
        result.append(orbit)
    return tuple(result)


def split_orbit_sensitivities(
    profile: Profile = (3, 1, 1, 5),
    step: float = 1e-5,
) -> tuple[SplitOrbitSensitivity, ...]:
    baseline_coefficients = forced_coefficients()
    baseline = fixed_beta_total(baseline_coefficients)
    result = []
    for orbit in split_orbits(profile):
        if all(entry in baseline_coefficients for entry in orbit):
            continue
        coefficients = dict(baseline_coefficients)
        for entry in orbit:
            coefficients.setdefault(entry, step)
        slope = (fixed_beta_total(coefficients) - baseline) / step
        if slope <= 1e-10:
            continue
        result.append(
            SplitOrbitSensitivity(
                entries=orbit,
                slope=slope,
                linear_threshold=(THRESHOLD - baseline) / slope,
            )
        )
    return tuple(sorted(result, key=lambda item: item.slope, reverse=True))


def exact_split_orbit_threshold(
    entries: tuple[ProfileSplit, ...],
) -> float:
    """Return the fixed-beta threshold when only one split orbit is added."""

    baseline_coefficients = forced_coefficients()

    def residual(value: float) -> float:
        coefficients = dict(baseline_coefficients)
        for entry in entries:
            coefficients.setdefault(entry, value)
        return fixed_beta_total(coefficients) - THRESHOLD

    upper = 1
    if residual(0) >= 0 or residual(upper) <= 0:
        raise AssertionError(
            ("split-orbit-threshold bracket", residual(0), residual(upper))
        )
    return float(brentq(residual, 0, upper, xtol=1e-13))


def reoptimized_split_orbit_threshold(
    entries: tuple[ProfileSplit, ...],
) -> tuple[float, float]:
    """Return the threshold after reoptimizing attenuation at each value."""

    def optimized_total(value: float) -> tuple[float, float]:
        return optimized_split_orbit_total(entries, value)

    root = brentq(
        lambda value: optimized_total(value)[0] - THRESHOLD,
        0,
        1,
        xtol=1e-12,
    )
    _, beta = optimized_total(root)
    return float(root), beta


def optimized_split_orbit_total(
    entries: tuple[ProfileSplit, ...],
    value: float,
) -> tuple[float, float]:
    """Reoptimize the repaired ledger after charging one split orbit."""

    coefficients = forced_coefficients()
    for entry in entries:
        coefficients.setdefault(entry, value)

    def total(beta: float) -> float:
        ledger = certificate(
            beta=beta,
            profile_split_coefficients=coefficients,
        )
        promise = euclidean_promise_concentration(DIMENSION, beta)
        return ledger.supporting_upper + promise.two_hypothesis_loss

    optimum = minimize_scalar(
        total,
        bounds=(0.75, 0.81),
        method="bounded",
        options={"xatol": 1e-12},
    )
    return float(optimum.fun), float(optimum.x)


def leading_scalar_obstruction(value: float) -> LeadingScalarObstruction:
    """Insert a proved lower witness for the leading split orbit."""

    total, beta = optimized_split_orbit_total(LEADING_SPLIT_ENTRIES, value)
    return LeadingScalarObstruction(
        coefficient=value,
        optimal_beta=beta,
        optimized_total=total,
        threshold_overshoot=total - THRESHOLD,
    )


def leading_hybrid_repair(value: float) -> LeadingHybridRepair:
    """Reoptimize the same physical cuts with hybrid promise packing."""

    coefficients = forced_coefficients()
    for entry in LEADING_SPLIT_ENTRIES:
        coefficients.setdefault(entry, value)

    def total(beta: float) -> float:
        ledger = certificate(
            beta=beta,
            profile_split_coefficients=coefficients,
        )
        promise = hybrid_euclidean_promise_concentration(DIMENSION, beta)
        return ledger.supporting_upper + promise.two_hypothesis_loss

    optimum = minimize_scalar(
        total,
        bounds=(0.75, 0.81),
        method="bounded",
        options={"xatol": 1e-12},
    )
    optimized = float(optimum.fun)
    return LeadingHybridRepair(
        coefficient=value,
        optimal_beta=float(optimum.x),
        optimized_total=optimized,
        threshold_slack=THRESHOLD - optimized,
    )


def hybrid_forced_coefficients() -> dict[ProfileSplit, float]:
    """Return both physical witness orbits at their certified coefficients."""

    result = forced_coefficients()
    for entry in LEADING_SPLIT_ENTRIES:
        result.setdefault(entry, CERTIFIED_LEADING_COEFFICIENT)
    return result


def fixed_hybrid_total(coefficients: dict[ProfileSplit, float]) -> float:
    ledger = certificate(
        beta=HYBRID_REPAIRED_BETA,
        profile_split_coefficients=coefficients,
    )
    promise = hybrid_euclidean_promise_concentration(
        DIMENSION,
        HYBRID_REPAIRED_BETA,
    )
    return ledger.supporting_upper + promise.two_hypothesis_loss


def hybrid_mapped_coefficients(
    profiles: tuple[Profile, ...],
    value: float,
) -> dict[ProfileSplit, float]:
    result = hybrid_forced_coefficients()
    for profile in profiles:
        for split in profile_splits(profile):
            result.setdefault((profile, split), value)
    return result


def hybrid_open_budget() -> RepairedOpenBudget:
    profiles = open_profiles()
    baseline = fixed_hybrid_total(hybrid_forced_coefficients())

    def residual(value: float) -> float:
        return (
            fixed_hybrid_total(hybrid_mapped_coefficients(profiles, value))
            - THRESHOLD
        )

    threshold = brentq(residual, 0, 1 / 10, xtol=1e-13)
    return RepairedOpenBudget(
        beta=HYBRID_REPAIRED_BETA,
        baseline_total=baseline,
        baseline_slack=THRESHOLD - baseline,
        common_coefficient_threshold=float(threshold),
    )


def extended_open_budget() -> RepairedOpenBudget:
    """Return the fixed-witness budget at the optimized two-split beta."""

    profiles = open_profiles()

    def total(coefficients: dict[ProfileSplit, float]) -> float:
        ledger = certificate(
            beta=EXTENDED_REPAIRED_BETA,
            profile_split_coefficients=coefficients,
        )
        promise = extended_euclidean_promise_concentration(
            DIMENSION,
            EXTENDED_REPAIRED_BETA,
        )
        return ledger.supporting_upper + promise.two_hypothesis_loss

    baseline = total(hybrid_forced_coefficients())

    def residual(value: float) -> float:
        return total(hybrid_mapped_coefficients(profiles, value)) - THRESHOLD

    threshold = brentq(residual, 0, 1 / 10, xtol=1e-13)
    return RepairedOpenBudget(
        beta=EXTENDED_REPAIRED_BETA,
        baseline_total=baseline,
        baseline_slack=THRESHOLD - baseline,
        common_coefficient_threshold=float(threshold),
    )


def hybrid_orbit_sensitivities(
    step: float = 1e-5,
) -> tuple[OrbitSensitivity, ...]:
    baseline = fixed_hybrid_total(hybrid_forced_coefficients())
    result = []
    for orbit in reversal_orbits(open_profiles()):
        slope = (
            fixed_hybrid_total(hybrid_mapped_coefficients(orbit, step))
            - baseline
        ) / step
        result.append(
            OrbitSensitivity(
                profiles=orbit,
                slope=slope,
                linear_threshold=(
                    (THRESHOLD - baseline) / slope
                    if slope > 1e-12
                    else float("inf")
                ),
            )
        )
    return tuple(sorted(result, key=lambda item: item.slope, reverse=True))


def hybrid_split_orbit_sensitivities(
    profile: Profile = (1, 1, 3, 5),
    step: float = 1e-5,
) -> tuple[SplitOrbitSensitivity, ...]:
    baseline_coefficients = hybrid_forced_coefficients()
    baseline = fixed_hybrid_total(baseline_coefficients)
    result = []
    for orbit in split_orbits(profile):
        if all(entry in baseline_coefficients for entry in orbit):
            continue
        coefficients = dict(baseline_coefficients)
        for entry in orbit:
            coefficients.setdefault(entry, step)
        slope = (fixed_hybrid_total(coefficients) - baseline) / step
        if slope <= 1e-10:
            continue
        result.append(
            SplitOrbitSensitivity(
                entries=orbit,
                slope=slope,
                linear_threshold=(THRESHOLD - baseline) / slope,
            )
        )
    return tuple(sorted(result, key=lambda item: item.slope, reverse=True))


def reoptimized_hybrid_split_threshold(
    entries: tuple[ProfileSplit, ...],
) -> tuple[float, float]:
    """Threshold one new split orbit under the hybrid promise theorem."""

    base = hybrid_forced_coefficients()

    def optimized_total(value: float) -> tuple[float, float]:
        coefficients = dict(base)
        for entry in entries:
            coefficients.setdefault(entry, value)

        def total(beta: float) -> float:
            ledger = certificate(
                beta=beta,
                profile_split_coefficients=coefficients,
            )
            promise = hybrid_euclidean_promise_concentration(DIMENSION, beta)
            return ledger.supporting_upper + promise.two_hypothesis_loss

        optimum = minimize_scalar(
            total,
            bounds=(0.75, 0.81),
            method="bounded",
            options={"xatol": 1e-12},
        )
        return float(optimum.fun), float(optimum.x)

    threshold = brentq(
        lambda value: optimized_total(value)[0] - THRESHOLD,
        0,
        1,
        xtol=1e-12,
    )
    _, beta = optimized_total(threshold)
    return float(threshold), beta


def optimized_extended_split_orbit_total(
    entries: tuple[ProfileSplit, ...],
    value: float,
) -> tuple[float, float]:
    """Reoptimize one split orbit with the two-split promise extension."""

    coefficients = hybrid_forced_coefficients()
    for entry in entries:
        coefficients.setdefault(entry, value)

    def total(beta: float) -> float:
        ledger = certificate(
            beta=beta,
            profile_split_coefficients=coefficients,
        )
        promise = extended_euclidean_promise_concentration(
            DIMENSION,
            beta,
        )
        return ledger.supporting_upper + promise.two_hypothesis_loss

    optimum = minimize_scalar(
        total,
        bounds=(0.75, 0.81),
        method="bounded",
        options={"xatol": 1e-12},
    )
    return float(optimum.fun), float(optimum.x)


def reoptimized_extended_split_threshold(
    entries: tuple[ProfileSplit, ...],
) -> tuple[float, float]:
    threshold = brentq(
        lambda value: (
            optimized_extended_split_orbit_total(entries, value)[0]
            - THRESHOLD
        ),
        0,
        1,
        xtol=1e-12,
    )
    _, beta = optimized_extended_split_orbit_total(entries, threshold)
    return float(threshold), beta


def extended_adjacent_diagnostic() -> ExtendedAdjacentDiagnostic:
    slices = horizontal_adjacent_slice_certificate(32)
    threshold, beta = reoptimized_extended_split_threshold(
        ADJACENT_SPLIT_ENTRIES
    )
    record_one_total, _ = optimized_extended_split_orbit_total(
        ADJACENT_SPLIT_ENTRIES,
        slices.record_one_coefficient,
    )
    combined_total, _ = optimized_extended_split_orbit_total(
        ADJACENT_SPLIT_ENTRIES,
        slices.combined_coefficient,
    )
    return ExtendedAdjacentDiagnostic(
        scalar_threshold=threshold,
        threshold_beta=beta,
        record_one_coefficient=slices.record_one_coefficient,
        record_three_coefficient=slices.record_three_coefficient,
        combined_slice_coefficient=slices.combined_coefficient,
        record_one_total=record_one_total,
        record_one_slack=THRESHOLD - record_one_total,
        combined_total=combined_total,
        combined_overshoot=combined_total - THRESHOLD,
    )


def adjacent_forced_coefficients() -> dict[ProfileSplit, float]:
    """Return all three certified physical witness orbits.

    The first two orbits are physical lower witnesses that must be charged
    at their certified values.  The adjacent orbit is charged at the safe
    direct record-one/record-three upper coefficient.  This last value is
    deliberately conservative; the accepted-sector repair makes a special
    compound improvement unnecessary for the current scalar target.
    """

    result = hybrid_forced_coefficients()
    adjacent = horizontal_adjacent_slice_certificate(32).combined_coefficient
    for entry in ADJACENT_SPLIT_ENTRIES:
        result.setdefault(entry, adjacent)
    return result


def coarse_open_completion_coefficients(
    value: float = 1 / 32,
) -> dict[ProfileSplit, float]:
    """Charge every other open high-degree split at one common value."""

    result = adjacent_forced_coefficients()
    for profile in open_profiles():
        for split in profile_splits(profile):
            result.setdefault((profile, split), value)
    return result


def optimized_extended_completion_total(value: float) -> tuple[float, float]:
    """Optimize attenuation for the coarse all-open completion target."""

    coefficients = coarse_open_completion_coefficients(value)

    def total(beta: float) -> float:
        ledger = certificate(
            beta=beta,
            profile_split_coefficients=coefficients,
        )
        promise = extended_euclidean_promise_concentration(DIMENSION, beta)
        return ledger.supporting_upper + promise.two_hypothesis_loss

    optimum = minimize_scalar(
        total,
        bounds=(0.75, 0.81),
        method="bounded",
        options={"xatol": 1e-12},
    )
    return float(optimum.fun), float(optimum.x)


def coarse_open_completion_target() -> CoarseOpenCompletionTarget:
    """Return the current scalar theorem target at q=32."""

    coefficient = 1 / 32
    forced = adjacent_forced_coefficients()
    coefficients = coarse_open_completion_coefficients(coefficient)
    total, beta = optimized_extended_completion_total(coefficient)
    return CoarseOpenCompletionTarget(
        coefficient=coefficient,
        optimal_beta=beta,
        optimized_total=total,
        threshold_slack=THRESHOLD - total,
        forced_entries=len(forced),
        coarse_entries=len(coefficients) - len(forced),
    )


def main() -> None:
    budget = repaired_open_budget()
    print(
        "repaired open-profile budget: "
        f"beta={budget.beta:.15g},"
        f"baseline_total={budget.baseline_total:.15g},"
        f"baseline_slack={budget.baseline_slack:.15g},"
        f"common_threshold={budget.common_coefficient_threshold:.15g},"
        f"inverse_threshold={1 / budget.common_coefficient_threshold:.12g}"
    )
    for item in orbit_sensitivities():
        print(
            "profile-orbit sensitivity: "
            f"profiles={item.profiles},"
            f"slope={item.slope:.12g},"
            f"linear_threshold={item.linear_threshold:.12g}"
        )
    split_ranking = split_orbit_sensitivities()
    for item in split_ranking:
        print(
            "leading split-orbit sensitivity: "
            f"entries={item.entries},"
            f"slope={item.slope:.12g},"
            f"linear_threshold={item.linear_threshold:.12g}"
        )
    print(
        "leading split-orbit exact fixed-beta threshold: "
        f"entries={split_ranking[0].entries},"
        f"threshold={exact_split_orbit_threshold(split_ranking[0].entries):.12g}"
    )
    threshold, beta = reoptimized_split_orbit_threshold(
        split_ranking[0].entries
    )
    print(
        "leading split-orbit reoptimized threshold: "
        f"threshold={threshold:.12g},beta={beta:.12g}"
    )
    hybrid = hybrid_open_budget()
    hybrid_profiles = hybrid_orbit_sensitivities()
    hybrid_splits = hybrid_split_orbit_sensitivities()
    hybrid_threshold, hybrid_beta = reoptimized_hybrid_split_threshold(
        hybrid_splits[0].entries
    )
    print(
        "hybrid repaired open-profile budget: "
        f"beta={hybrid.beta:.15g},"
        f"baseline_total={hybrid.baseline_total:.15g},"
        f"baseline_slack={hybrid.baseline_slack:.15g},"
        f"common_threshold={hybrid.common_coefficient_threshold:.15g},"
        f"top_profile={hybrid_profiles[0].profiles},"
        f"top_split={hybrid_splits[0].entries},"
        f"top_split_threshold={hybrid_threshold:.12g},"
        f"top_split_beta={hybrid_beta:.12g}"
    )
    adjacent = extended_adjacent_diagnostic()
    extended = extended_open_budget()
    print(
        "extended adjacent-split diagnostic: "
        f"baseline_beta={extended.beta:.12g},"
        f"baseline_total={extended.baseline_total:.12g},"
        f"baseline_slack={extended.baseline_slack:.12g},"
        f"common_threshold={extended.common_coefficient_threshold:.12g},"
        f"threshold={adjacent.scalar_threshold:.12g},"
        f"threshold_beta={adjacent.threshold_beta:.12g},"
        f"record_one={adjacent.record_one_coefficient:.12g},"
        f"record_three={adjacent.record_three_coefficient:.12g},"
        f"combined={adjacent.combined_slice_coefficient:.12g},"
        f"record_one_slack={adjacent.record_one_slack:.12g},"
        f"combined_overshoot={adjacent.combined_overshoot:.12g}"
    )
    coarse = coarse_open_completion_target()
    print(
        "coarse all-open completion target: "
        f"coefficient={coarse.coefficient:.12g},"
        f"beta={coarse.optimal_beta:.12g},"
        f"total={coarse.optimized_total:.12g},"
        f"slack={coarse.threshold_slack:.12g},"
        f"forced_entries={coarse.forced_entries},"
        f"coarse_entries={coarse.coarse_entries}"
    )


if __name__ == "__main__":
    main()
