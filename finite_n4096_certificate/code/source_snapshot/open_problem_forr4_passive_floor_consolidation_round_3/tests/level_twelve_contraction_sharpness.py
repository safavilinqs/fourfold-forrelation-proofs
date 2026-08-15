#!/usr/bin/env python3
"""Regression for the level-twelve distinctness-masked sharpness witness."""

from __future__ import annotations

from fractions import Fraction
from math import isclose
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from level_twelve_contraction_sharpness import (  # noqa: E402
    direct_masked_contraction,
    grouped_test_vector_norms,
    lower_witness_ratio_to_n_inverse,
    masked_lower_witness_exact,
    sharpness_audit,
    valid_positive_history_count,
)


def main() -> None:
    expected_values = {
        4: Fraction(9, 512),
        8: Fraction(2646, 65536),
        16: Fraction(617400, 16777216),
    }
    for dimension, expected in expected_values.items():
        exact = masked_lower_witness_exact(dimension)
        if exact != expected:
            raise AssertionError(("exact masked value", dimension, exact, expected))
        direct = direct_masked_contraction(dimension)
        if not isclose(direct, float(exact), rel_tol=0.0, abs_tol=3e-13):
            raise AssertionError(("direct masked value", dimension, direct, exact))
        norms = grouped_test_vector_norms(dimension)
        if len(norms) != 11 or any(
            not isclose(norm, 1.0, rel_tol=0.0, abs_tol=3e-13) for norm in norms
        ):
            raise AssertionError(("unit grouped vectors", dimension, norms))

    if valid_positive_history_count() != 180:
        raise AssertionError("level-twelve positive history count")
    audit = sharpness_audit(1024)
    if (
        audit.initial_configurations != 12
        or audit.terminal_histories != 1080
        or audit.displayed_history_orders != 180
        or not audit.positive_coefficient_family
    ):
        raise AssertionError(("positive terminal family", audit))
    if audit.legal_assigned_sigma != 1:
        raise AssertionError(("legal paired placement", audit))
    if audit.projective_upper_exponent != -1 or audit.masked_lower_exponent != -1:
        raise AssertionError(("matching upper/lower exponents", audit))
    if audit.masked_lower_value != masked_lower_witness_exact(1024):
        raise AssertionError(("N=1024 lower value", audit))
    if audit.ratio_to_n_inverse != lower_witness_ratio_to_n_inverse(1024):
        raise AssertionError(("N=1024 N^-1 ratio", audit))
    if audit.ratio_to_n_inverse <= Fraction(99, 100):
        raise AssertionError(("distinctness should cost below one percent", audit))
    if -audit.projective_upper_exponent / 12 != Fraction(1, 12):
        raise AssertionError(("current exponent barrier", audit))

    print(
        "level-twelve contraction sharpness passed: "
        f"initials={audit.initial_configurations},"
        f"histories={audit.terminal_histories},"
        f"displayed_orders={audit.displayed_history_orders},"
        f"positive={audit.positive_coefficient_family},"
        f"assigned_sigma={audit.legal_assigned_sigma},"
        f"upper_exponent={audit.projective_upper_exponent},"
        f"lower_exponent={audit.masked_lower_exponent},"
        f"N1024_lower={float(audit.masked_lower_value):.12g},"
        f"N1024_ratio={float(audit.ratio_to_n_inverse):.12g},"
        f"barrier={-audit.projective_upper_exponent / 12}"
    )


if __name__ == "__main__":
    main()
