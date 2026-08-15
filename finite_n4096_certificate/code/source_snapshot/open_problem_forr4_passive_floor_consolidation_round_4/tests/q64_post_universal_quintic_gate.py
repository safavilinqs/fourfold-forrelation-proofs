#!/usr/bin/env python3
"""Regression for the live post-universal q64 quintic gate."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_post_universal_quintic_gate import (  # noqa: E402
    artifact_text,
    diagnostic,
    quintic_entries,
    quintic_split_depth,
)


def main() -> None:
    result = diagnostic()
    observed_discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.proved_entries,
        result.remaining_open_entries,
        result.quintic_entries,
        result.extreme_split_entries,
        result.balanced_split_entries,
    )
    expected_discrete = (64, 4096, 16_384, 220, 668, 184, 104, 80)
    if observed_discrete != expected_discrete:
        raise AssertionError(("post-universal quintic discrete result", observed_discrete))

    observed = (
        result.current_routing.total,
        result.current_routing.margin_to_one_third,
        result.common_reserve_coefficient,
        result.extreme_local_slice_coefficient,
        result.balanced_local_slice_coefficient,
        result.balanced_reserve_coefficient_after_extreme_local_slice,
        result.local_slice_proxy.total,
        result.local_slice_proxy.margin_to_one_third,
        result.local_slice_proxy_overshoot,
    )
    expected = (
        0.33193582943438404,
        0.0013975038989492705,
        0.12526109565056054,
        0.12397463639031407,
        0.14955611574342903,
        0.12568133975135265,
        0.33824808166480846,
        -0.0049147483314751494,
        0.0049147483314751494,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-10):
        raise AssertionError(("post-universal quintic numeric result", observed))

    depths = [quintic_split_depth(entry) for entry in quintic_entries()]
    if depths.count(1) != 104 or depths.count(2) != 80:
        raise AssertionError(("quintic split partition", depths))
    if not (
        result.extreme_local_slice_coefficient < result.common_reserve_coefficient
        < result.balanced_local_slice_coefficient
    ):
        raise AssertionError("quintic local-scale decision ordering")

    committed = (
        ROOT / "artifacts" / "q64_post_universal_quintic_gate.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 post-universal quintic artifact")

    print(
        "q64 post-universal quintic gate passed: "
        f"entries={result.quintic_entries},"
        f"common_reserve={result.common_reserve_coefficient:.12g},"
        "balanced_reserve="
        f"{result.balanced_reserve_coefficient_after_extreme_local_slice:.12g},"
        f"proxy_overshoot={result.local_slice_proxy_overshoot:.12g},"
        "status=routing_gate_not_theorem"
    )


if __name__ == "__main__":
    main()
