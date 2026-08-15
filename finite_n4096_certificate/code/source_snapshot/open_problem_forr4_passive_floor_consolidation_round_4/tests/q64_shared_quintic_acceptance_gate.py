#!/usr/bin/env python3
"""Regression for the q64 shared-quintic and adaptive acceptance gate."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_shared_quintic_acceptance_gate import (  # noqa: E402
    artifact_text,
    diagnostic,
    remaining_quintic_entries,
)


def main() -> None:
    result = diagnostic()
    discrete = (
        result.order,
        result.dimension,
        result.proved_entries,
        result.remaining_open_entries,
        result.remaining_quintic_entries,
        result.extreme_quintic_entries,
        result.balanced_quintic_entries,
        result.remaining_nonquintic_entries,
        result.higher_split_only_in_cubic_profile_entries,
        result.noncubic_profile_entries,
        result.two_split_cubics_one_split_higher_entries,
        result.one_split_cubic_no_split_higher_entries,
    )
    if discrete != (64, 4096, 380, 508, 48, 32, 16, 460, 176, 140, 96, 48):
        raise AssertionError(("shared-quintic discrete gate", discrete))
    observed = (
        result.current_routing.total,
        result.common_reserve_coefficient,
        result.extreme_sufficient_coefficient,
        result.balanced_sufficient_coefficient,
        result.balanced_reserve_gate_after_extreme,
        result.sufficient_two_tier_proxy.total,
        result.sufficient_two_tier_proxy.margin_to_one_third,
        result.raw_adaptive_overhead_cap,
        result.adaptive_overhead_cap_retaining_allowance,
        result.raw_adaptive_multiplier_cap,
        result.adaptive_multiplier_cap_retaining_allowance,
    )
    expected = (
        0.32605080644648143,
        0.41031455336741013,
        0.12397463639031407,
        0.14955611574342903,
        1.0933655828876938,
        0.3262031888677698,
        0.007130144465563537,
        0.007130144465563537,
        0.006130144465563537,
        1.0218579851727134,
        1.0187924112171953,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("shared-quintic numeric gate", observed))
    if len(remaining_quintic_entries()) != 48:
        raise AssertionError("shared-quintic entry inventory")
    committed = (
        ROOT / "artifacts" / "q64_shared_quintic_acceptance_gate.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 shared-quintic acceptance artifact")
    print(
        "q64 shared quintic acceptance gate passed: "
        f"entries={result.remaining_quintic_entries},"
        f"common_reserve={result.common_reserve_coefficient:.12g},"
        f"two_tier_total={result.sufficient_two_tier_proxy.total:.12g},"
        "adaptive_additive_cap="
        f"{result.adaptive_overhead_cap_retaining_allowance:.12g},"
        "status=routing_and_adaptive_requirement_not_theorem"
    )


if __name__ == "__main__":
    main()
