#!/usr/bin/env python3
"""Regression for the exact alternating double-endpoint Walsh certificate."""

from __future__ import annotations

from decimal import Decimal
import importlib.util
from pathlib import Path
import sys

import numpy as np


SEARCHES = Path(__file__).resolve().parents[1] / "searches"
sys.path.insert(0, str(SEARCHES))

from alternating_double_endpoint_spectrum import (
    pair_arrays,
    scaled_endpoint_weights,
    spectrum_certificate,
)


def load_direct_benchmark():
    path = (
        Path(__file__).resolve().parents[2]
        / "open_problem_forr4_passive_floor_consolidation_round_2"
        / "searches"
        / "double_endpoint_joint_schur_benchmark.py"
    )
    spec = importlib.util.spec_from_file_location(
        "round_two_double_endpoint_benchmark", path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(("cannot load direct benchmark", path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_endpoint_formula() -> None:
    benchmark = load_direct_benchmark()
    for order in (2, 4):
        endpoint, _ = benchmark.endpoint_slice(order)
        dimension = order * order
        hadamard = (
            np.kron(
                benchmark.sylvester(order),
                benchmark.sylvester(order),
            ).astype(float)
            / order
        )
        pair_left, pair_right = pair_arrays(dimension)
        pair_xor = pair_left ^ pair_right
        formula = np.empty_like(endpoint)
        for singleton in range(dimension):
            weights = (
                scaled_endpoint_weights(
                    order, singleton, pair_left, pair_right
                ).astype(float)
                / (order - 1)
            )
            formula[singleton] = (
                hadamard[singleton ^ pair_xor].T * weights
            )
        if not np.array_equal(endpoint, formula):
            raise AssertionError(("endpoint moment formula", order))


def main() -> None:
    check_endpoint_formula()
    expected = {
        2: Decimal("0.47159181589114324"),
        4: Decimal("0.06420087162467479"),
    }
    for order, direct_value in expected.items():
        certificate = spectrum_certificate(order)
        if abs(certificate.coefficient - direct_value) > Decimal("3e-16"):
            raise AssertionError(
                ("direct row-Gram mismatch", order, certificate.coefficient)
            )
        if int(certificate.multiplicities.sum()) != (order * order) ** 2:
            raise AssertionError(("block rank", order))

    q32 = spectrum_certificate(32)
    if q32.denominator != 968381956096:
        raise AssertionError(("q32 denominator", q32.denominator))
    if len(q32.numerators) != 20:
        raise AssertionError(("q32 spectrum classes", len(q32.numerators)))
    if int(q32.numerators[0]) != 889351684:
        raise AssertionError(("q32 minimum numerator", q32.numerators[0]))
    if int(q32.numerators[-1]) != 7945127554564:
        raise AssertionError(("q32 maximum numerator", q32.numerators[-1]))
    target = Decimal(
        "0.00035759351713982535053208591764918627335980824852312"
    )
    if abs(q32.coefficient - target) > Decimal("3e-49"):
        raise AssertionError(("q32 coefficient", q32.coefficient))

    print(
        "alternating double-endpoint spectrum passed: "
        f"q2={expected[2]},q4={expected[4]},"
        f"q32={q32.coefficient},classes={len(q32.numerators)},"
        f"rank={int(q32.multiplicities.sum())}"
    )


if __name__ == "__main__":
    main()
