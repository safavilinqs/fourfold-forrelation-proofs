#!/usr/bin/env python3
"""Regression for the coefficient-one distinctness-mask audit."""

from __future__ import annotations

from pathlib import Path
import sys
from collections import Counter

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_masked_universal_audit import (  # noqa: E402
    Q2_SPLITS,
    affected_classes,
    affected_entries,
    artifact_text,
    diagnostic,
    exact_q2_matrix,
    theorem_registry,
)


def main() -> None:
    result = diagnostic()
    classes = affected_classes()
    if tuple(len(entries) for entries in classes) != (96, 14, 24, 124, 96):
        raise AssertionError("coefficient-one class inventory")
    if len(affected_entries()) != 354:
        raise AssertionError("coefficient-one distinct inventory")
    discrete = (
        result.affected_entries,
        result.masked_quintic_slice_repaired_entries,
        result.masked_local_walsh_repaired_entries,
        result.masked_cubic_endpoint_repaired_entries,
        result.masked_double_quintic_endpoint_repaired_entries,
        result.masked_double_quintic_record_repaired_entries,
        result.masked_four_cubic_incidence_repaired_entries,
        result.masked_cubic_septimic_chain_repaired_entries,
        result.masked_recovered_cubic_quintic_incidence_repaired_entries,
        result.masked_joint_recovered_cubic_quintic_repaired_entries,
        result.final_residual_chain_repaired_entries,
        result.total_masked_repaired_entries,
        result.remaining_quarantined_entries,
        result.supported_entries_before_dual_caveat,
        result.open_entries_before_dual_caveat,
        result.dual_endpoint_certified_entries,
        result.dual_endpoint_caveat_entries,
        result.conservative_supported_entries,
        result.conservative_open_entries,
        result.abstract_mask_counterexample_dimension,
    )
    if discrete != (
        354,
        54,
        180,
        12,
        6,
        12,
        38,
        12,
        28,
        12,
        80,
        354,
        0,
        888,
        0,
        12,
        0,
        888,
        0,
        4,
    ):
        raise AssertionError(("masked-universal registry", discrete))
    statuses = Counter(row["status"] for row in theorem_registry())
    if statuses != {
        "proved_nonuniversal_inherited": 442,
        "proved_dual_endpoint_schur": 12,
        "proved_masked_quintic_slice": 54,
        "proved_masked_local_walsh": 180,
        "proved_masked_cubic_endpoint": 12,
        "proved_masked_double_quintic_endpoint": 6,
        "proved_masked_double_quintic_record": 12,
        "proved_masked_four_cubic_incidence": 38,
        "proved_masked_cubic_septimic_chain": 12,
        "proved_masked_recovered_cubic_quintic_incidence": 28,
        "proved_masked_joint_recovered_cubic_quintic": 12,
        "proved_final_residual_chain": 80,
    }:
        raise AssertionError(("masked-universal registry statuses", statuses))
    if not np.isclose(result.abstract_mask_counterexample_value, 1.5):
        raise AssertionError("J-I coefficient counterexample")

    expected_shapes = (
        (96, 96, 2304),
        (256, 216, 6912),
        (576, 576, 20736),
    )
    expected_uniform = (
        0.41666666666666663,
        0.5192148880023929,
        0.6106334218568518,
    )
    expected_upper = (
        0.416666666766667,
        0.5843669705203611,
        0.6804748216924796,
    )
    for index, (split, row) in enumerate(zip(Q2_SPLITS, result.q2_rows)):
        integer = exact_q2_matrix(split)
        if set(np.unique(integer)) != {-1, 0, 1}:
            raise AssertionError(("q2 exact entry alphabet", split))
        if (row.rows, row.columns, row.nonzero_entries) != expected_shapes[index]:
            raise AssertionError(("q2 matrix shape", split))
        if row.scaled_entry_denominator != 8:
            raise AssertionError(("q2 moment denominator", split))
        if not np.isclose(
            row.uniform_weighted_nuclear,
            expected_uniform[index],
            rtol=0,
            atol=3e-13,
        ):
            raise AssertionError(("q2 uniform nuclear", split))
        if not np.isclose(
            row.uniform_supporting_upper,
            expected_upper[index],
            rtol=0,
            atol=3e-12,
        ):
            raise AssertionError(("q2 tangent upper", split))
        if row.uniform_supporting_upper >= 1:
            raise AssertionError(("q2 coefficient-one survival", split))

    observed = (
        result.complete_floating_routing.total,
        result.complete_floating_routing.margin_to_one_third,
        result.generic_mask_minimum_coefficient,
        result.generic_mask_maximum_coefficient,
        result.generic_mask_repaired_routing.total,
        result.generic_mask_repaired_routing.margin_to_one_third,
    )
    expected = (
        0.23991166777023,
        0.0934216655631033,
        3.0,
        33.970562748477136,
        3.4034411220503005,
        -3.070107788716967,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-11):
        raise AssertionError(("masked-universal q64 routing", observed))

    committed = (
        ROOT / "artifacts" / "q64_masked_universal_audit.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale masked-universal audit artifact")
    print(
        "q64 masked-universal audit passed: "
        f"affected={result.affected_entries},"
        f"supported={result.supported_entries_before_dual_caveat},"
        f"conservative_supported={result.conservative_supported_entries},"
        "q2_upper="
        + "/".join(
            f"{row.uniform_supporting_upper:.12g}" for row in result.q2_rows
        )
        + f",generic_mask_total={result.generic_mask_repaired_routing.total:.12g},"
        "status=all_888_registry_entries_certified_one_batch"
    )


if __name__ == "__main__":
    main()
