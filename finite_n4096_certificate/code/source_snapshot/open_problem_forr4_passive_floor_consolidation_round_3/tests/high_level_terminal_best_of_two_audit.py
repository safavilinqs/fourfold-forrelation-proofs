#!/usr/bin/env python3
"""Regression for the exact high-level terminal best-of-two audit."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from high_level_terminal_best_of_two_audit import (  # noqa: E402
    enumerate_terminal_states,
    high_level_terminal_audit,
    sigma_one_audit,
)


EXPECTED_FAILURES = (
    (
        9,
        (9,),
        (((0, 0), (0, 1)), ((0, 0), (0, 1), (1, 2)), ((0, 0), (1, 1), (2, 2))),
    ),
    (
        9,
        (9,),
        (((0, 0), (1, 1), (2, 2)), ((0, 0), (1, 0), (2, 1)), ((0, 0), (1, 0))),
    ),
    (
        10,
        (6, 4),
        (((0, 0), (1, 1), (2, 2)), ((0, 0), (1, 0), (2, 1)), ((0, 0), (1, 1))),
    ),
)


def main() -> None:
    audit = high_level_terminal_audit()
    if audit.reachable_states != 222 or audit.terminal_states != 34:
        raise AssertionError(("exact state counts", audit))
    if audit.high_level_terminals != 17 or audit.sensitive_high_level_terminals != 8:
        raise AssertionError(("high-level counts", audit))
    if audit.terminals_by_level != ((9, 8), (10, 6), (11, 2), (12, 1)):
        raise AssertionError(("level counts", audit.terminals_by_level))
    if audit.sensitive_by_level != ((9, 3), (10, 3), (11, 1), (12, 1)):
        raise AssertionError(("sensitive level counts", audit.sensitive_by_level))
    if (
        audit.sigma_one_types,
        audit.passing_sigma_one_types,
        audit.failing_sigma_one_types,
    ) != (
        8,
        5,
        3,
    ):
        raise AssertionError(("best-of-two verdict counts", audit))
    if audit.worst_best_decay != Fraction(
        1, 2
    ) or audit.worst_required_decay != Fraction(5, 8):
        raise AssertionError(("worst high-level score", audit))
    if audit.worst_safe_decay_by_level != (
        (9, Fraction(1, 2)),
        (10, Fraction(1, 2)),
        (11, Fraction(1)),
        (12, Fraction(1)),
    ):
        raise AssertionError(("levelwise safe decay", audit))
    if audit.previous_global_exponent != Fraction(1, 24):
        raise AssertionError(("previous global exponent", audit))
    if audit.proved_global_exponent != Fraction(1, 20):
        raise AssertionError(("improved global exponent", audit))

    observed_failures = tuple(
        (failure.level, failure.component_sizes, failure.boundaries)
        for failure in audit.failures
    )
    if observed_failures != EXPECTED_FAILURES:
        raise AssertionError(("joint saturator list", observed_failures))
    for failure in audit.failures:
        if failure.projective_decay != Fraction(1, 2):
            raise AssertionError(("joint saturator projective decay", failure))
        if failure.best_decay != Fraction(1, 2) or failure.passes:
            raise AssertionError(("joint saturator verdict", failure))

    _, terminals = enumerate_terminal_states()
    high = [terminal for terminal in terminals if terminal.level >= 9]
    for terminal in high:
        for component in terminal.components:
            if {layer for layer, _ in component} != {0, 1, 2, 3}:
                raise AssertionError(("high-level component misses a layer", terminal))

    level_twelve = [terminal for terminal in high if terminal.level == 12]
    if len(level_twelve) != 1:
        raise AssertionError(("unique level-twelve type", level_twelve))
    terminal = level_twelve[0]
    if not terminal.sensitive or tuple(map(len, terminal.components)) != (4, 4, 4):
        raise AssertionError(("level-twelve three-path forest", terminal))
    terminal_audit = sigma_one_audit(terminal)
    if terminal_audit is None:
        raise AssertionError("level-twelve forest should admit assigned sigma one")
    if (
        terminal_audit.projective_decay != 1
        or terminal_audit.required_decay != Fraction(3, 4)
    ):
        raise AssertionError(("level-twelve projective repair", terminal_audit))

    print(
        "high-level terminal best-of-two audit passed: "
        f"reachable={audit.reachable_states},"
        f"terminals={audit.terminal_states},"
        f"high={audit.high_level_terminals},"
        f"sensitive={audit.sensitive_high_level_terminals},"
        f"sigma_one={audit.sigma_one_types},"
        f"passing={audit.passing_sigma_one_types},"
        f"joint_saturators={audit.failing_sigma_one_types},"
        f"failure_levels={tuple(failure.level for failure in audit.failures)},"
        f"level12_projective={terminal_audit.projective_decay},"
        f"global_exponent={audit.proved_global_exponent}"
    )


if __name__ == "__main__":
    main()
