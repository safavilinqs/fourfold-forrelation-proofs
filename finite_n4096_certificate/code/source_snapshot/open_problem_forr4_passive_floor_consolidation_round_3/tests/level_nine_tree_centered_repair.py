#!/usr/bin/env python3
"""Regression checks for the centered repair of both level-nine trees."""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from level_nine_tree_centered_repair import (  # noqa: E402
    A0,
    B0,
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
    if history.new_transfers != (3,) or history.existing_transfers != (2,):
        raise AssertionError(("transfer types", history))
    if history.derivative_events != (2,):
        raise AssertionError(("derivative events", history))
    if history.derivative_sites != (("A", A0), ("A", B0)):
        raise AssertionError(("centered derivative sites", history))
    expected_weights = (
        (("A", A0), (LocalWeight("gamma", 1),)),
        (
            ("A", B0),
            (LocalWeight("gamma", 1), LocalWeight("stein", 1)),
        ),
    )
    if history.derivative_weights_by_site != expected_weights:
        raise AssertionError(("centered derivative weights", history))
    if history.time_exponents != (1, 2, 2):
        raise AssertionError(("time exponents", history))

    repair = repair_audit()
    if repair.dangerous_partitions != 7:
        raise AssertionError(("dangerous partitions", repair))
    type_a_counts = (
        repair.marked_neighbor_partitions,
        repair.existing_middle_partitions,
        repair.existing_endpoint_partitions,
        repair.fresh_endpoint_partitions,
    )
    if type_a_counts != (14, 132, 855, 1905):
        raise AssertionError(("Type-A branch counts", repair))
    type_a_decays = (
        repair.minimum_marked_neighbor_decay,
        repair.minimum_existing_middle_decay,
        repair.minimum_fresh_endpoint_decay,
    )
    if type_a_decays != (Fraction(1),) * 3:
        raise AssertionError(("Type-A unit decays", repair))
    if repair.minimum_existing_endpoint_decay != Fraction(3, 2):
        raise AssertionError(("Type-A endpoint decay", repair))

    if repair.reflected_dangerous_partitions != 7:
        raise AssertionError(("reflected dangerous partitions", repair))
    type_b_counts = (
        repair.reflected_marked_neighbor_partitions,
        repair.reflected_existing_middle_partitions,
        repair.reflected_existing_endpoint_partitions,
        repair.reflected_cancelled_fresh_partitions,
    )
    if type_b_counts != (14, 132, 855, 1905):
        raise AssertionError(("Type-B branch counts", repair))
    type_b_decays = (
        repair.reflected_minimum_marked_neighbor_decay,
        repair.reflected_minimum_existing_middle_decay,
        repair.reflected_minimum_existing_endpoint_decay,
    )
    if type_b_decays != (Fraction(1), Fraction(1), Fraction(3, 2)):
        raise AssertionError(("Type-B retained decays", repair))
    if repair.proved_global_exponent != Fraction(1, 16):
        raise AssertionError(("global exponent", repair))

    print(
        "level-nine reflected-tree centered repair passed: "
        f"histories={history.histories},"
        f"profiles={history.weight_profiles},"
        f"derivative_sites={history.derivative_sites},"
        f"type_a_counts={type_a_counts},"
        f"type_b_counts={type_b_counts},"
        f"global_exponent={repair.proved_global_exponent}"
    )


if __name__ == "__main__":
    main()
