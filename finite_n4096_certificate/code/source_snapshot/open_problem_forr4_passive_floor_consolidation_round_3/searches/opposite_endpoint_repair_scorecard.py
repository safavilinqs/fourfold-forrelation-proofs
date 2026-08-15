#!/usr/bin/env python3
"""Rank scalar-ledger repairs after the mixed-orbit obstruction.

For each accepted profile family, uniformly scale all of its proved
fixed-split coefficients, reoptimize attenuation, and solve for the smallest
reduction that offsets the forced opposite-endpoint cuts.  This is a route
selection diagnostic, not a claim that the reductions have been proved.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log

from scipy.optimize import brentq, minimize_scalar

from attenuation_promise_concentration import promise_concentration
from occupation_compatible_sector_optimization import (
    SEPARATED_CUBIC_QUINTIC_PROFILES,
    SEPARATED_QUINTIC_CUBIC_PROFILES,
    TRIPLE_CUBIC_PROFILES,
    certificate,
    coefficient,
    profile_splits,
)
from opposite_endpoint_vertical_mixture_witness import (
    forced_split_coefficients,
    vertical_mixture_witness,
)


@dataclass(frozen=True)
class RepairTarget:
    name: str
    profiles: tuple[tuple[int, ...], ...]
    required_reduction: float
    optimal_beta: float


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


def optimized_total(
    coefficients: dict[tuple[tuple[int, ...], tuple[int, ...]], float],
) -> tuple[float, float]:
    def total(beta: float) -> float:
        ledger = certificate(
            beta=beta,
            profile_split_coefficients=coefficients,
        )
        return (
            ledger.supporting_upper
            + promise_concentration(1024, beta).two_hypothesis_loss
        )

    result = minimize_scalar(
        total,
        bounds=(0.77, 0.79),
        method="bounded",
        options={"xatol": 2e-9},
    )
    return float(result.fun), float(result.x)


def repair_target(
    name: str,
    profiles: tuple[tuple[int, ...], ...],
    forced: dict[tuple[tuple[int, ...], tuple[int, ...]], float],
) -> RepairTarget:
    proved = {
        (profile, split): coefficient(profile, split, True, None)
        for profile in profiles
        for split in profile_splits(profile)
    }

    def coefficients(scale: float):
        result = dict(forced)
        result.update(
            {key: scale * value for key, value in proved.items()}
        )
        return result

    if optimized_total(coefficients(0))[0] >= 1 / 3:
        raise AssertionError(("profile family cannot repair budget", name))
    scale = brentq(
        lambda value: optimized_total(coefficients(value))[0] - 1 / 3,
        0,
        1,
        xtol=2e-7,
    )
    _, beta = optimized_total(coefficients(scale))
    return RepairTarget(
        name=name,
        profiles=profiles,
        required_reduction=1 - scale,
        optimal_beta=beta,
    )


def promise_proxy_target(beta: float, overshoot: float) -> tuple[float, float]:
    promise = promise_concentration(1024, beta)
    target_loss = promise.two_hypothesis_loss - overshoot
    target_proxy = -promise.promise_gap**2 / (2 * log(target_loss / 2))
    return target_proxy, 1 - target_proxy / promise.four_chain_proxy


def main() -> None:
    witness = vertical_mixture_witness(32)
    forced = forced_split_coefficients(witness.coefficient)
    forced_total, forced_beta = optimized_total(forced)
    overshoot = forced_total - 1 / 3
    proxy, proxy_reduction = promise_proxy_target(forced_beta, overshoot)
    print(
        "promise repair target: "
        f"beta={forced_beta:.12g},"
        f"overshoot={overshoot:.12g},"
        f"target_proxy={proxy:.12g},"
        f"proxy_reduction={proxy_reduction:.12g}"
    )
    targets = sorted(
        (
            repair_target(name, profiles, forced)
            for name, profiles in GROUPS.items()
        ),
        key=lambda result: result.required_reduction,
    )
    for target in targets:
        print(
            "profile repair target: "
            f"name={target.name},"
            f"required_reduction={target.required_reduction:.12g},"
            f"optimal_beta={target.optimal_beta:.12g}"
        )


if __name__ == "__main__":
    main()
