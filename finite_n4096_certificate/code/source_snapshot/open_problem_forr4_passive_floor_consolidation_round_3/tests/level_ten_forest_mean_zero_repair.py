#!/usr/bin/env python3
"""Regression checks for the centered repair of the level-ten forest."""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from level_ten_forest_mean_zero_repair import (  # noqa: E402
    C0,
    LocalWeight,
    enumerate_histories,
    repair_audit,
)


def main() -> None:
    history = enumerate_histories()
    if history.potential_twelve_initial_configurations != 6:
        raise AssertionError(("initial configurations", history))
    if history.contributing_initial_configurations != 4:
        raise AssertionError(("contributing initial configurations", history))
    if history.histories != 200 or history.weight_profiles != 4:
        raise AssertionError(("history enumeration", history))
    if history.new_transfers != (4,) or history.existing_transfers != (1,):
        raise AssertionError(("transfer types", history))
    if history.derivative_events != (1,):
        raise AssertionError(("derivative events", history))
    if history.common_derivative_site != ("B", C0):
        raise AssertionError(("centered derivative site", history))
    expected_weights = (
        LocalWeight("gamma", 1),
        LocalWeight("stein", 1),
    )
    if history.common_derivative_weights != expected_weights:
        raise AssertionError(("centered derivative weights", history))
    if history.time_exponents != (2, 2, 1):
        raise AssertionError(("time exponents", history))

    repair = repair_audit()
    if repair.dangerous_partitions != 282:
        raise AssertionError(("dangerous partitions", repair))
    if repair.duplicate_branches != 564 or repair.bridge_partitions != 282:
        raise AssertionError(("marked-neighbor branches", repair))
    if repair.new_neighbor_partitions != 5295:
        raise AssertionError(("new-neighbor partitions", repair))
    decays = (
        repair.minimum_duplicate_decay,
        repair.minimum_bridge_decay,
        repair.minimum_existing_outer_decay,
        repair.repaired_level_ten_decay,
        repair.repaired_level_eleven_decay,
    )
    if decays != (Fraction(1),) * 5:
        raise AssertionError(("repaired decay", repair))
    if repair.proved_global_exponent != Fraction(1, 18):
        raise AssertionError(("global exponent", repair))

    print(
        "level-ten forest mean-zero repair passed: "
        f"potential_initials={history.potential_twelve_initial_configurations},"
        f"contributing_initials={history.contributing_initial_configurations},"
        f"histories={history.histories},"
        f"profiles={history.weight_profiles},"
        f"centered_weights={history.common_derivative_weights},"
        f"partitions={repair.dangerous_partitions},"
        f"extended={repair.new_neighbor_partitions},"
        f"duplicate_decay={repair.minimum_duplicate_decay},"
        f"bridge_decay={repair.minimum_bridge_decay},"
        f"outer_decay={repair.minimum_existing_outer_decay},"
        f"global_exponent={repair.proved_global_exponent}"
    )


if __name__ == "__main__":
    main()
