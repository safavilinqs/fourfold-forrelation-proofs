#!/usr/bin/env python3
"""Regression for the exact projective translation-cocycle inventory."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from masked_translation_cocycle_inventory import (  # noqa: E402
    artifact_text,
    diagnostic,
    predicted_normalized_rank,
)


EXPECTED_RANKS = {
    4: {0: 21, 8: 26, 16: 50},
    8: {0: 21, 12: 26, 24: 50},
    64: {0: 21, 24: 26, 48: 50},
}
EXPECTED_NORMALIZED = {0: 21, 4: 26, 8: 50}
EXPECTED_PARITIES = {
    "0000:0": 15,
    "0011:4": 19,
    "0101:8": 19,
    "0110:8": 19,
    "1001:8": 5,
    "1010:8": 7,
    "1100:4": 7,
    "1111:0": 6,
}


def main() -> None:
    for order in (4, 8, 64):
        result = diagnostic(order)
        if result.templates != 97:
            raise AssertionError(("template count", order, result.templates))
        if result.rank_counts != EXPECTED_RANKS[order]:
            raise AssertionError(("commutator ranks", order, result.rank_counts))
        if result.normalized_rank_counts != EXPECTED_NORMALIZED:
            raise AssertionError(("normalized ranks", order))
        if result.parity_rank_counts != EXPECTED_PARITIES:
            raise AssertionError(("parity classification", order))
        for row in result.rows:
            if row.row_commutator_rank != row.column_commutator_rank:
                raise AssertionError(("row/column rank mismatch", order, row))
            if not row.row_column_cocycles_equal:
                raise AssertionError(("row/column cocycle mismatch", order, row))
            if row.normalized_rank != predicted_normalized_rank(
                row.split_parities
            ):
                raise AssertionError(("closed rank formula", order, row))
        committed = (
            ROOT
            / "artifacts"
            / f"q{order}_masked_translation_cocycle_inventory.json"
        ).read_text(encoding="utf-8")
        if committed != artifact_text(result):
            raise AssertionError(("stale cocycle inventory artifact", order))
    print(
        "masked translation-cocycle inventories passed: "
        "q4_ranks=0/8/16,q8_ranks=0/12/24,q64_ranks=0/24/48,"
        "normalized_types=0/4/8,status=exact_three_type_classification"
    )


if __name__ == "__main__":
    main()
