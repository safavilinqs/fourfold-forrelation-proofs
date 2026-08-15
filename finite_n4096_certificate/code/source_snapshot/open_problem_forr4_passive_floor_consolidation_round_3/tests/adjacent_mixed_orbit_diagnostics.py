#!/usr/bin/env python3
"""Exact small-q checks for the adjacent mixed-orbit diagnostics."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from adjacent_cubic_quintic_mixed_orbit_q4 import (
    build_data,
    direct_evaluate_small_support,
    evaluate,
    pure_laws,
)
from adjacent_vertical_triple_symmetric_witness import (
    closed_symmetric_vertical_triple_witness,
    symmetric_vertical_triple_witness,
)


def main() -> None:
    data = build_data(4)
    row, column = pure_laws(data, 4, 1, (0, 4, 1))
    pure = evaluate(data, row, column).objective
    if not np.isclose(pure, 17 / 576, atol=2e-14):
        raise AssertionError(("pure adjacent orbit", pure))

    # A genuinely mixed law checks the Fourier reduction against an
    # independently assembled occurrence matrix.
    row = np.zeros_like(row)
    row[0, 3] = 0.4
    row[4, 1] = 0.6
    column = np.zeros_like(column)
    column[0] = 0.55
    column[1] = 0.45
    fourier = evaluate(data, row, column).objective
    direct = direct_evaluate_small_support(data, row, column)
    if not np.isclose(fourier, direct, atol=3e-13):
        raise AssertionError(("mixed adjacent reduction", fourier, direct))

    # The 25-frequency-class closed formula must reproduce the independent
    # complete q=4 construction before it is used at q=32.
    exact_symmetric = symmetric_vertical_triple_witness(4)
    closed_symmetric = closed_symmetric_vertical_triple_witness(4)
    if not np.isclose(
        exact_symmetric.coefficient,
        closed_symmetric.coefficient,
        atol=3e-13,
    ):
        raise AssertionError(
            ("closed vertical-triple formula", exact_symmetric, closed_symmetric)
        )
    if not np.isclose(
        closed_symmetric.coefficient,
        0.185024058902,
        atol=3e-12,
    ):
        raise AssertionError(("vertical-triple q4 witness", closed_symmetric))

    print(
        "adjacent mixed-orbit diagnostics passed: "
        f"pure={pure:.12g},"
        f"mixed={fourier:.12g},"
        f"mixed_direct={direct:.12g},"
        f"symmetric_q4={closed_symmetric.coefficient:.12g}"
    )


if __name__ == "__main__":
    main()
