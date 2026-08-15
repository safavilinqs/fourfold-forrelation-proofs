#!/usr/bin/env python3
"""Regression for the complete dependency-exact outward q64 ledger."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from json import loads
from math import exp
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_complete_outward_ledger import (  # noqa: E402
    BETA,
    RESERVE_THRESHOLD,
    artifact_text,
    coefficient_upper_map,
    complete_coefficients,
    diagnostic,
    promise_upper,
    registry_rows,
)
from q64_dual_endpoint_schur_insertion import (  # noqa: E402
    cubic_fixed_pair_energy,
    dual_endpoint_entries,
    quintic_fixed_triple_energy,
)
from q64_final_residual_chain_contraction import (  # noqa: E402
    repaired_entries as final_entries,
    squared_coefficient as final_squared_coefficient,
)
from q64_paper_target_gate import balanced_open_entries  # noqa: E402


def main() -> None:
    result = diagnostic()
    raw = complete_coefficients()
    upper = coefficient_upper_map()
    rows = registry_rows()
    if BETA != Fraction(19, 25):
        raise AssertionError("ledger beta changed")
    if len(raw) != 888 or len(rows) != 888:
        raise AssertionError(("high-sector dependency registry", len(raw), len(rows)))
    if len(balanced_open_entries()) != 888:
        raise AssertionError("balanced inventory")
    if any(2 * sum(split) != sum(profile) for profile, split in raw):
        raise AssertionError("unbalanced routing placeholder entered theorem ledger")
    if any(upper[entry] <= Fraction.from_float(value) for entry, value in raw.items()):
        raise AssertionError("coefficient grid is not strictly outward")

    dual_squared = cubic_fixed_pair_energy() * quintic_fixed_triple_energy()
    if any(upper[entry] ** 2 < dual_squared for entry in dual_endpoint_entries()):
        raise AssertionError("dual exact coefficient is not dominated")
    if any(
        upper[entry] ** 2 < final_squared_coefficient(entry)
        for entry in final_entries()
    ):
        raise AssertionError("final residual exact coefficient is not dominated")

    statuses = Counter(row["status"] for row in rows)
    expected = {
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
    }
    if statuses != expected:
        raise AssertionError(("balanced dependency statuses", statuses))

    exponent, promise = promise_upper()
    if float(promise) < 2 * exp(-float(exponent)):
        raise AssertionError("rational promise relaxation rounded inward")
    total = Fraction(str(result.total_upper))
    threshold = Fraction(str(result.reserve_threshold))
    if threshold > RESERVE_THRESHOLD or total >= RESERVE_THRESHOLD:
        raise AssertionError(("reserve gate", total, threshold))
    if not result.passes_reserve_gate:
        raise AssertionError("reserve decision")
    if not (Fraction(1, 4) < total < Fraction(27, 100)):
        raise AssertionError(("unexpected complete total", total))
    if result.supported_balanced_entries != 888 or result.open_balanced_entries:
        raise AssertionError("complete registry count")
    if (
        result.high_sector_profile_splits,
        result.certified_balanced_high_sector_coefficients,
        result.excluded_unbalanced_high_sector_entries,
        result.excluded_unbalanced_high_sector_incidence_records,
        result.excluded_unbalanced_high_sector_undirected_edges,
    ) != (6016, 888, 5128, 272, 136):
        raise AssertionError("number-sector scope audit")

    committed_path = ROOT / "artifacts" / "q64_complete_outward_ledger.json"
    committed = committed_path.read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale complete outward ledger artifact")
    payload = loads(committed)
    if len(payload["coefficient_registry"]) != 888:
        raise AssertionError("artifact dependency registry")
    print(
        "q64 complete outward ledger passed: "
        f"certified_balanced_coefficients="
        f"{result.certified_balanced_high_sector_coefficients},"
        "excluded_unbalanced_incidences=272/136_edges,"
        f"supported={result.supported_balanced_entries},"
        f"perron={result.collatz_perron_upper},"
        f"promise={result.promise_loss_upper},"
        f"total={result.total_upper},"
        f"reserve_margin={result.reserve_margin_lower},"
        "status=eligible_for_adaptive_lift"
    )


if __name__ == "__main__":
    main()
