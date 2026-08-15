#!/usr/bin/env python3
"""Regression for the universal q64 septimic theorem insertion."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_block_coherent_contraction import block_coherent_entries  # noqa: E402
from q64_chain_aware_insertion import chain_aware_entries  # noqa: E402
from q64_remaining_class_gates import contraction_class  # noqa: E402
from q64_universal_septimic_insertion import (  # noqa: E402
    artifact_text,
    diagnostic,
    septimic_entries,
)


def main() -> None:
    result = diagnostic()
    observed_discrete = (
        result.order,
        result.dimension,
        result.sign_modes,
        result.septimic_entries,
        result.previous_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
    )
    expected_discrete = (64, 4096, 16_384, 96, 110, 206, 682)
    if observed_discrete != expected_discrete:
        raise AssertionError(("universal septimic discrete result", observed_discrete))

    observed = (
        result.universal_coefficient,
        result.previous_routing.total,
        result.septimic_inserted.total,
        result.septimic_inserted.beta,
        result.septimic_inserted.perron_upper,
        result.septimic_inserted.promise_loss,
        result.septimic_inserted.margin_to_one_third,
        result.margin_spent,
        result.reserve_after_declared_allowance,
    )
    expected = (
        1.0,
        0.2960908671821436,
        0.3293832216221608,
        0.7461075872914213,
        0.3121610995953498,
        0.017222122026810997,
        0.003950111711172488,
        0.03329235444001721,
        0.002950111711172488,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("universal septimic numeric result", observed))

    entries = septimic_entries()
    if len(entries) != 96 or not all(
        7 in profile
        and contraction_class(entry) == "one_split_cubic_one_split_higher"
        for entry in entries
        for profile, _ in (entry,)
    ):
        raise AssertionError("universal septimic entry class")
    overlap = set(entries).intersection(block_coherent_entries())
    overlap.update(set(entries).intersection(chain_aware_entries()))
    if overlap:
        raise AssertionError(("septimic theorem overlap", overlap))

    committed = (
        ROOT / "artifacts" / "q64_universal_septimic_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale q64 universal septimic artifact")

    print(
        "q64 universal septimic insertion passed: "
        f"entries={result.septimic_entries},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.septimic_inserted.total:.12g},"
        f"margin={result.septimic_inserted.margin_to_one_third:.12g},"
        "status=quarantined_unmasked_coefficient_one_diagnostic"
    )


if __name__ == "__main__":
    main()
