#!/usr/bin/env python3
"""Regression for the exact-order q64 finite-size routing gate."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_paper_target_gate import artifact_text, diagnostic  # noqa: E402
import degree_six_joint_occupation_optimization as degree_six  # noqa: E402
import double_endpoint_occupation_optimization as double_endpoint  # noqa: E402
import occupation_compatible_sector_optimization as occupation  # noqa: E402


def main() -> None:
    inherited_orders = (
        occupation.ORDER,
        degree_six.ORDER,
        double_endpoint.ORDER,
    )
    result = diagnostic()
    restored_orders = (
        occupation.ORDER,
        degree_six.ORDER,
        double_endpoint.ORDER,
    )
    if restored_orders != inherited_orders or inherited_orders != (32, 32, 32):
        raise AssertionError(
            ("order context restoration", inherited_orders, restored_orders)
        )

    observed_discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.balanced_open_entries,
        result.cubic_profile_entries,
        result.noncubic_profile_entries,
        result.block_coherent_entries,
        result.internally_split_cubic_entries,
        result.cubic_profile_high_internal_entries,
        result.noncubic_high_internal_entries,
    )
    if observed_discrete != (
        64,
        4096,
        16_384,
        888,
        724,
        164,
        70,
        486,
        192,
        140,
    ):
        raise AssertionError(("q64 discrete target", observed_discrete))

    observed = (
        result.common_one_over_q.total,
        result.common_one_over_q.beta,
        result.common_one_over_sqrt_q.total,
        result.common_one_over_sqrt_q.beta,
        result.common_threshold_coefficient,
        result.common_reserve_coefficient,
        result.cubic_slice_target,
        result.noncubic_target,
        result.two_tier_target.total,
        result.two_tier_target.beta,
        result.two_tier_target.perron_upper,
        result.two_tier_target.promise_loss,
        result.two_tier_target.margin_to_one_third,
    )
    expected = (
        0.11349030834996607,
        0.7520370567856766,
        0.2418944198500511,
        0.7477797317355201,
        0.1999106655420832,
        0.1990891760722882,
        0.12403521525363623,
        0.5,
        0.3191811621612196,
        0.7463284992438618,
        0.3028473373687459,
        0.016333824792473674,
        0.014152171172113703,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("q64 numeric target", observed))
    if result.common_reserve_coefficient <= 1 / np.sqrt(64):
        raise AssertionError("q64 reserve gate no longer covers the square-root target")
    if result.two_tier_target.margin_to_one_third <= 0.01:
        raise AssertionError("q64 two-tier routing reserve became too small")

    committed = (ROOT / "artifacts" / "q64_paper_target_gate.json").read_text(
        encoding="utf-8"
    )
    if committed != artifact_text(result):
        raise AssertionError("stale q64 paper-target artifact")

    print(
        "q64 paper target gate passed: "
        f"N={result.dimension},"
        f"M={result.sign_modes},"
        f"common_gate={result.common_threshold_coefficient:.12g},"
        f"reserve_gate={result.common_reserve_coefficient:.12g},"
        f"two_tier_total={result.two_tier_target.total:.12g},"
        f"two_tier_margin={result.two_tier_target.margin_to_one_third:.12g},"
        "status=routing_target_not_theorem"
    )


if __name__ == "__main__":
    main()
