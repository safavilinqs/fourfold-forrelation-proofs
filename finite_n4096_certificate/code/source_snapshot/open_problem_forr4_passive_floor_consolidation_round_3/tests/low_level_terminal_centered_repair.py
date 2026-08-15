#!/usr/bin/env python3
"""Regression checks for the full low-level terminal centered repair."""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from level_ten_forest_mean_zero_repair import LocalWeight  # noqa: E402
from low_level_terminal_centered_repair import complete_audit  # noqa: E402


def main() -> None:
    audit = complete_audit()
    image = audit.image
    if image.initial_potentials != (12, 8, 8, 4):
        raise AssertionError(("initial potentials", image))
    if image.reachable_states != 236:
        raise AssertionError(("reachable states", image))
    if image.terminal_types != 39 or image.sensitive_terminal_types != 22:
        raise AssertionError(("terminal counts", image))
    expected_terminal_counts = (
        (4, 1),
        (5, 2),
        (6, 4),
        (7, 6),
        (8, 9),
        (9, 8),
        (10, 6),
        (11, 2),
        (12, 1),
    )
    expected_sensitive_counts = (
        (4, 1),
        (5, 2),
        (6, 3),
        (7, 4),
        (8, 4),
        (9, 3),
        (10, 3),
        (11, 1),
        (12, 1),
    )
    if image.terminal_types_by_level != expected_terminal_counts:
        raise AssertionError(("terminal levels", image))
    if image.sensitive_types_by_level != expected_sensitive_counts:
        raise AssertionError(("sensitive levels", image))
    expected_decays = (
        (4, Fraction(1, 2)),
        (5, Fraction(1)),
        (6, Fraction(1, 2)),
        (7, Fraction(1, 2)),
        (8, Fraction(1)),
        (9, Fraction(1, 2)),
        (10, Fraction(1, 2)),
        (11, Fraction(1)),
        (12, Fraction(1)),
    )
    if image.minimum_safe_decay_by_level != expected_decays:
        raise AssertionError(("safe decay levels", image))
    if image.sensitive_level_eight_types != 4:
        raise AssertionError(("level-eight types", image))
    if image.level_eight_minimum_decay != Fraction(1):
        raise AssertionError(("level-eight decay", image))
    if image.level_seven_saturators != 1 or image.level_six_saturators != 1:
        raise AssertionError(("low-level saturators", image))

    seven = audit.level_seven_coefficients
    if seven.potential_initial_configurations != ((8, 4), (12, 2)):
        raise AssertionError(("level-seven initials", seven))
    if seven.contributing_initial_configurations != ((8, 4), (12, 2)):
        raise AssertionError(("level-seven contributing initials", seven))
    if seven.histories_by_initial_potential != ((8, 12), (12, 24)):
        raise AssertionError(("level-seven history potentials", seven))
    if seven.histories != 36 or seven.weight_profiles != 8:
        raise AssertionError(("level-seven histories", seven))
    if seven.fresh_transfers_by_initial_potential != ((8, (2,)), (12, (1,))):
        raise AssertionError(("level-seven fresh transfers", seven))
    if seven.existing_transfers_by_initial_potential != (
        (8, (1,)),
        (12, (2,)),
    ):
        raise AssertionError(("level-seven existing transfers", seven))
    if seven.derivative_events_by_initial_potential != (
        (8, (1,)),
        (12, (1,)),
    ):
        raise AssertionError(("level-seven derivatives", seven))
    if seven.common_derivative_sites != (("A", (0, 0)),):
        raise AssertionError(("level-seven centered site", seven))
    if seven.common_derivative_weights != (
        (("A", (0, 0)), (LocalWeight("gamma", 1),)),
    ):
        raise AssertionError(("level-seven centered weight", seven))
    if seven.time_exponents != (1, 1, 1):
        raise AssertionError(("level-seven time exponents", seven))

    six = audit.level_six_coefficients
    if six.potential_initial_configurations != ((8, 2),):
        raise AssertionError(("level-six initials", six))
    if six.contributing_initial_configurations != ((8, 2),):
        raise AssertionError(("level-six contributing initials", six))
    if six.histories_by_initial_potential != ((8, 4),):
        raise AssertionError(("level-six history potentials", six))
    if six.histories != 4 or six.weight_profiles != 2:
        raise AssertionError(("level-six histories", six))
    if six.fresh_transfers_by_initial_potential != ((8, (1,)),):
        raise AssertionError(("level-six fresh transfers", six))
    if six.existing_transfers_by_initial_potential != ((8, (1,)),):
        raise AssertionError(("level-six existing transfers", six))
    if six.derivative_events_by_initial_potential != ((8, (1,)),):
        raise AssertionError(("level-six derivatives", six))
    if six.common_derivative_sites != (("A", (1, 0)),):
        raise AssertionError(("level-six centered site", six))
    if six.common_derivative_weights != ((("A", (1, 0)), (LocalWeight("gamma", 1),)),):
        raise AssertionError(("level-six centered weight", six))
    if six.time_exponents != (0, 1, 1):
        raise AssertionError(("level-six time exponents", six))

    repair = audit.repair
    level_seven_counts = (
        repair.level_seven_marked_neighbor_cases,
        repair.level_seven_existing_middle_cases,
        repair.level_seven_existing_endpoint_cases,
        repair.level_seven_fresh_endpoint_cases,
    )
    if repair.level_seven_dangerous_partitions != 10:
        raise AssertionError(("level-seven partitions", repair))
    if level_seven_counts != (20, 106, 586, 1692):
        raise AssertionError(("level-seven branch counts", repair))
    level_seven_decays = (
        repair.level_seven_minimum_marked_neighbor_decay,
        repair.level_seven_minimum_existing_middle_decay,
        repair.level_seven_minimum_existing_endpoint_decay,
        repair.level_seven_minimum_fresh_endpoint_decay,
    )
    if level_seven_decays != (
        Fraction(1),
        Fraction(1),
        Fraction(3, 2),
        Fraction(1),
    ):
        raise AssertionError(("level-seven decay", repair))
    if repair.level_six_dangerous_partitions != 31:
        raise AssertionError(("level-six partitions", repair))
    if repair.level_six_marked_outer_cases != 31:
        raise AssertionError(("level-six marked branches", repair))
    if repair.level_six_cancelled_fresh_outer_cases != 144:
        raise AssertionError(("level-six cancelled branches", repair))
    if repair.level_six_minimum_marked_outer_decay != Fraction(1):
        raise AssertionError(("level-six marked decay", repair))
    if repair.proved_global_exponent != Fraction(1, 12):
        raise AssertionError(("global exponent", repair))

    print(
        "low-level terminal centered repair passed: "
        f"reachable={image.reachable_states},"
        f"terminals={image.terminal_types},"
        f"sensitive={image.sensitive_terminal_types},"
        f"level8_decay={image.level_eight_minimum_decay},"
        f"level7_histories={seven.histories},"
        f"level7_counts={level_seven_counts},"
        f"level6_histories={six.histories},"
        f"level6_cancelled={repair.level_six_cancelled_fresh_outer_cases},"
        f"global_exponent={repair.proved_global_exponent}"
    )


if __name__ == "__main__":
    main()
