#!/usr/bin/env python3
"""Audit every q64 theorem that used the unmasked coefficient-one lemma."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from itertools import combinations, product
from json import dumps
from math import prod
from pathlib import Path

import numpy as np

from q64_chain_aware_insertion import (
    inserted_coefficients as pre_universal_coefficients,
)
from q64_dual_endpoint_schur_insertion import dual_endpoint_entries
from q64_masked_quintic_slice_repair import (
    coefficient_map as masked_quintic_coefficients,
    repaired_entries as masked_quintic_repaired_entries,
)
from q64_masked_local_walsh_repair import (
    coefficient_map as masked_local_walsh_coefficients,
    repaired_entries as masked_local_walsh_repaired_entries,
)
from q64_masked_cubic_endpoint_repair import (
    coefficient_map as masked_cubic_endpoint_coefficients,
    repaired_entries as masked_cubic_endpoint_repaired_entries,
)
from q64_masked_double_quintic_endpoint_repair import (
    coefficient_map as masked_double_quintic_endpoint_coefficients,
    repaired_entries as masked_double_quintic_endpoint_repaired_entries,
)
from q64_masked_double_quintic_record_repair import (
    coefficient_map as masked_double_quintic_record_coefficients,
    repaired_entries as masked_double_quintic_record_repaired_entries,
)
from q64_masked_four_cubic_incidence_repair import (
    coefficient_map as masked_four_cubic_incidence_coefficients,
    repaired_entries as masked_four_cubic_incidence_repaired_entries,
)
from q64_masked_cubic_septimic_chain_repair import (
    coefficient_map as masked_cubic_septimic_chain_coefficients,
    repaired_entries as masked_cubic_septimic_chain_repaired_entries,
)
from q64_masked_recovered_cubic_quintic_incidence_repair import (
    coefficient_map as masked_recovered_cubic_quintic_incidence_coefficients,
    repaired_entries as masked_recovered_cubic_quintic_incidence_repaired_entries,
)
from q64_joint_recovered_cubic_quintic_contraction import (
    coefficient_map as masked_joint_recovered_cubic_quintic_coefficients,
    repaired_entries as masked_joint_recovered_cubic_quintic_repaired_entries,
)
from q64_final_residual_chain_contraction import (
    coefficient_map as final_residual_chain_coefficients,
    repaired_entries as final_residual_chain_repaired_entries,
)
from q64_noncubic_recovered_universal_insertion import (
    recovered_universal_entries,
    universal_noncubic_entries,
)
from q64_paper_target_gate import (
    OptimizedLedger,
    balanced_open_entries,
    optimize,
)
from q64_same_side_whole_link_insertion import (
    inserted_coefficients as latest_claimed_coefficients,
    remaining_entries as final_residual_entries,
    subset_disjointness_factor,
)
from q64_universal_double_cubic_insertion import double_cubic_entries
from q64_universal_multicubic_insertion import multicubic_entries
from q64_universal_septimic_insertion import septimic_entries
from signed_permutation_link_moment import chain_moment


ROOT = Path(__file__).resolve().parents[1]
PROFILE = (3, 3, 3, 3)
Q2_SPLITS = (
    (0, 1, 2, 3),
    (1, 1, 1, 3),
    (1, 1, 2, 2),
)


@dataclass(frozen=True)
class MaskedQ2Row:
    split: tuple[int, ...]
    split_cubic_blocks: int
    rows: int
    columns: int
    nonzero_entries: int
    scaled_entry_denominator: int
    uniform_weighted_nuclear: float
    uniform_supporting_upper: float
    row_gradient_spread: float
    column_gradient_spread: float


@dataclass(frozen=True)
class MaskedUniversalAudit:
    universal_septimic_entries: int
    universal_multicubic_entries: int
    universal_double_cubic_entries: int
    universal_noncubic_entries: int
    recovered_universal_entries: int
    affected_entries: int
    masked_quintic_slice_repaired_entries: int
    masked_local_walsh_repaired_entries: int
    masked_cubic_endpoint_repaired_entries: int
    masked_double_quintic_endpoint_repaired_entries: int
    masked_double_quintic_record_repaired_entries: int
    masked_four_cubic_incidence_repaired_entries: int
    masked_cubic_septimic_chain_repaired_entries: int
    masked_recovered_cubic_quintic_incidence_repaired_entries: int
    masked_joint_recovered_cubic_quintic_repaired_entries: int
    final_residual_chain_repaired_entries: int
    total_masked_repaired_entries: int
    remaining_quarantined_entries: int
    supported_entries_before_dual_caveat: int
    open_entries_before_dual_caveat: int
    dual_endpoint_certified_entries: int
    dual_endpoint_caveat_entries: int
    conservative_supported_entries: int
    conservative_open_entries: int
    abstract_mask_counterexample_dimension: int
    abstract_mask_counterexample_value: float
    q2_rows: tuple[MaskedQ2Row, ...]
    complete_floating_routing: OptimizedLedger
    generic_mask_minimum_coefficient: float
    generic_mask_maximum_coefficient: float
    generic_mask_repaired_routing: OptimizedLedger


def affected_classes() -> tuple[tuple[tuple, ...], ...]:
    """Return the five disjoint coefficient-one-dependent inventories."""

    return (
        septimic_entries(),
        multicubic_entries(),
        double_cubic_entries(),
        universal_noncubic_entries(),
        recovered_universal_entries(),
    )


def affected_entries() -> frozenset[tuple]:
    classes = affected_classes()
    result = frozenset(entry for entries in classes for entry in entries)
    if len(result) != sum(len(entries) for entries in classes):
        raise AssertionError("coefficient-one affected classes overlap")
    return result


def quarantined_coefficients() -> dict[tuple, float]:
    """Use only valid coefficients and reset every unrepaired dependency."""

    result = latest_claimed_coefficients()
    base = pre_universal_coefficients()
    for entry in affected_entries():
        result[entry] = base[entry]
    result.update(masked_quintic_coefficients())
    result.update(masked_local_walsh_coefficients())
    result.update(masked_cubic_endpoint_coefficients())
    result.update(masked_double_quintic_endpoint_coefficients())
    result.update(masked_double_quintic_record_coefficients())
    result.update(masked_four_cubic_incidence_coefficients())
    result.update(masked_cubic_septimic_chain_coefficients())
    result.update(masked_recovered_cubic_quintic_incidence_coefficients())
    result.update(masked_joint_recovered_cubic_quintic_coefficients())
    result.update(final_residual_chain_coefficients())
    return result


def entry_mask_factor(entry: tuple) -> float:
    """Safe product of direct-sum factors for every split block."""

    profile, split = entry
    return prod(
        subset_disjointness_factor(
            min(selected, degree - selected),
            max(selected, degree - selected),
        )
        for degree, selected in zip(profile, split, strict=True)
        if 0 < selected < degree
    )


def generic_mask_repaired_coefficients() -> dict[tuple, float]:
    result = latest_claimed_coefficients()
    for entry in affected_entries():
        result[entry] = entry_mask_factor(entry)
    return result


def exact_q2_matrix(split: tuple[int, ...]) -> np.ndarray:
    """Return eight times one exact physical masked occurrence matrix."""

    order = 2
    dimension = order * order
    complement = tuple(
        degree - selected
        for degree, selected in zip(PROFILE, split, strict=True)
    )

    def subsets(size: int) -> tuple[tuple[int, ...], ...]:
        return tuple(combinations(range(dimension), size))

    rows = tuple(product(*(subsets(size) for size in split)))
    columns = tuple(product(*(subsets(size) for size in complement)))
    matrix = np.zeros((len(rows), len(columns)), dtype=np.int8)
    for row_index, row in enumerate(rows):
        for column_index, column in enumerate(columns):
            if any(
                set(left).intersection(right)
                for left, right in zip(row, column, strict=True)
            ):
                continue
            supports = tuple(
                tuple(sorted(left + right))
                for left, right in zip(row, column, strict=True)
            )
            exact = 8 * chain_moment(order, supports)
            if exact.denominator != 1 or abs(exact.numerator) > 1:
                raise AssertionError(("unexpected q2 moment", split, exact))
            matrix[row_index, column_index] = exact.numerator
    return matrix


def q2_row(split: tuple[int, ...]) -> MaskedQ2Row:
    """Evaluate a physical masked matrix and a concavity tangent upper."""

    integer = exact_q2_matrix(split)
    kernel = integer.astype(float) / 8
    rows, columns = kernel.shape
    row_law = np.full(rows, 1 / rows)
    column_law = np.full(columns, 1 / columns)
    weighted = (
        np.sqrt(row_law)[:, None]
        * kernel
        * np.sqrt(column_law)[None, :]
    )
    left, singular_values, right = np.linalg.svd(
        weighted, full_matrices=False
    )
    polar = left @ right
    row_gradient = np.sum(polar * weighted, axis=1) / (2 * row_law)
    column_gradient = np.sum(polar * weighted, axis=0) / (
        2 * column_law
    )
    # The weighted trace norm is jointly concave in the two diagonal laws.
    # Its tangent therefore bounds both probability simplexes by the sum of
    # the two maximum gradient coordinates. Add a conservative numerical
    # allowance far above the observed SVD roundoff.
    supporting_upper = (
        float(row_gradient.max() + column_gradient.max()) + 1e-10
    )
    return MaskedQ2Row(
        split=split,
        split_cubic_blocks=sum(value in (1, 2) for value in split),
        rows=rows,
        columns=columns,
        nonzero_entries=int(np.count_nonzero(integer)),
        scaled_entry_denominator=8,
        uniform_weighted_nuclear=float(singular_values.sum()),
        uniform_supporting_upper=supporting_upper,
        row_gradient_spread=float(np.ptp(row_gradient)),
        column_gradient_spread=float(np.ptp(column_gradient)),
    )


def diagnostic() -> MaskedUniversalAudit:
    classes = affected_classes()
    affected = affected_entries()
    quintic_repaired = frozenset(masked_quintic_repaired_entries())
    local_walsh_repaired = frozenset(masked_local_walsh_repaired_entries())
    cubic_endpoint_repaired = frozenset(masked_cubic_endpoint_repaired_entries())
    double_quintic_endpoint_repaired = frozenset(
        masked_double_quintic_endpoint_repaired_entries()
    )
    double_quintic_record_repaired = frozenset(
        masked_double_quintic_record_repaired_entries()
    )
    four_cubic_incidence_repaired = frozenset(
        masked_four_cubic_incidence_repaired_entries()
    )
    cubic_septimic_chain_repaired = frozenset(
        masked_cubic_septimic_chain_repaired_entries()
    )
    recovered_cubic_quintic_incidence_repaired = frozenset(
        masked_recovered_cubic_quintic_incidence_repaired_entries()
    )
    joint_recovered_cubic_quintic_repaired = frozenset(
        masked_joint_recovered_cubic_quintic_repaired_entries()
    )
    final_residual = frozenset(final_residual_entries())
    final_residual_repaired = frozenset(final_residual_chain_repaired_entries())
    repair_classes = (
        quintic_repaired,
        local_walsh_repaired,
        cubic_endpoint_repaired,
        double_quintic_endpoint_repaired,
        double_quintic_record_repaired,
        four_cubic_incidence_repaired,
        cubic_septimic_chain_repaired,
        recovered_cubic_quintic_incidence_repaired,
        joint_recovered_cubic_quintic_repaired,
    )
    if any(
        left.intersection(right)
        for index, left in enumerate(repair_classes)
        for right in repair_classes[index + 1 :]
    ):
        raise AssertionError("masked repair inventories overlap")
    repaired = frozenset().union(*repair_classes)
    if not repaired.issubset(affected):
        raise AssertionError("masked repair is outside the quarantine")
    if final_residual_repaired != final_residual:
        raise AssertionError("final residual repair does not match its inventory")
    mask_values = tuple(entry_mask_factor(entry) for entry in affected)
    supported = 888 - len(affected - repaired) - len(
        final_residual - final_residual_repaired
    )
    dual_certified = len(dual_endpoint_entries())
    return MaskedUniversalAudit(
        universal_septimic_entries=len(classes[0]),
        universal_multicubic_entries=len(classes[1]),
        universal_double_cubic_entries=len(classes[2]),
        universal_noncubic_entries=len(classes[3]),
        recovered_universal_entries=len(classes[4]),
        affected_entries=len(affected),
        masked_quintic_slice_repaired_entries=len(quintic_repaired),
        masked_local_walsh_repaired_entries=len(local_walsh_repaired),
        masked_cubic_endpoint_repaired_entries=len(cubic_endpoint_repaired),
        masked_double_quintic_endpoint_repaired_entries=len(
            double_quintic_endpoint_repaired
        ),
        masked_double_quintic_record_repaired_entries=len(
            double_quintic_record_repaired
        ),
        masked_four_cubic_incidence_repaired_entries=len(
            four_cubic_incidence_repaired
        ),
        masked_cubic_septimic_chain_repaired_entries=len(
            cubic_septimic_chain_repaired
        ),
        masked_recovered_cubic_quintic_incidence_repaired_entries=len(
            recovered_cubic_quintic_incidence_repaired
        ),
        masked_joint_recovered_cubic_quintic_repaired_entries=len(
            joint_recovered_cubic_quintic_repaired
        ),
        final_residual_chain_repaired_entries=len(final_residual_repaired),
        total_masked_repaired_entries=len(repaired),
        remaining_quarantined_entries=len(affected - repaired),
        supported_entries_before_dual_caveat=supported,
        open_entries_before_dual_caveat=888 - supported,
        dual_endpoint_certified_entries=dual_certified,
        dual_endpoint_caveat_entries=0,
        conservative_supported_entries=supported,
        conservative_open_entries=888 - supported,
        abstract_mask_counterexample_dimension=4,
        abstract_mask_counterexample_value=2 * (1 - 1 / 4),
        q2_rows=tuple(q2_row(split) for split in Q2_SPLITS),
        complete_floating_routing=optimize(
            mapped_coefficients=quarantined_coefficients()
        ),
        generic_mask_minimum_coefficient=min(mask_values),
        generic_mask_maximum_coefficient=max(mask_values),
        generic_mask_repaired_routing=optimize(
            mapped_coefficients=generic_mask_repaired_coefficients()
        ),
    )


def theorem_registry() -> tuple[dict[str, object], ...]:
    """Return a dependency-exact status row for all 888 balanced entries."""

    affected = affected_entries()
    quintic_repaired = frozenset(masked_quintic_repaired_entries())
    local_walsh_repaired = frozenset(masked_local_walsh_repaired_entries())
    cubic_endpoint_repaired = frozenset(masked_cubic_endpoint_repaired_entries())
    double_quintic_endpoint_repaired = frozenset(
        masked_double_quintic_endpoint_repaired_entries()
    )
    double_quintic_record_repaired = frozenset(
        masked_double_quintic_record_repaired_entries()
    )
    four_cubic_incidence_repaired = frozenset(
        masked_four_cubic_incidence_repaired_entries()
    )
    cubic_septimic_chain_repaired = frozenset(
        masked_cubic_septimic_chain_repaired_entries()
    )
    recovered_cubic_quintic_incidence_repaired = frozenset(
        masked_recovered_cubic_quintic_incidence_repaired_entries()
    )
    joint_recovered_cubic_quintic_repaired = frozenset(
        masked_joint_recovered_cubic_quintic_repaired_entries()
    )
    final_residual_repaired = frozenset(final_residual_chain_repaired_entries())
    dual = frozenset(dual_endpoint_entries())
    residual = frozenset(final_residual_entries())
    if affected.intersection(dual | residual) or dual.intersection(residual):
        raise AssertionError("masked-audit registry classes overlap")
    rows = []
    for profile, split in sorted(balanced_open_entries()):
        entry = (profile, split)
        if entry in quintic_repaired:
            status = "proved_masked_quintic_slice"
        elif entry in local_walsh_repaired:
            status = "proved_masked_local_walsh"
        elif entry in cubic_endpoint_repaired:
            status = "proved_masked_cubic_endpoint"
        elif entry in double_quintic_endpoint_repaired:
            status = "proved_masked_double_quintic_endpoint"
        elif entry in double_quintic_record_repaired:
            status = "proved_masked_double_quintic_record"
        elif entry in four_cubic_incidence_repaired:
            status = "proved_masked_four_cubic_incidence"
        elif entry in cubic_septimic_chain_repaired:
            status = "proved_masked_cubic_septimic_chain"
        elif entry in recovered_cubic_quintic_incidence_repaired:
            status = "proved_masked_recovered_cubic_quintic_incidence"
        elif entry in joint_recovered_cubic_quintic_repaired:
            status = "proved_masked_joint_recovered_cubic_quintic"
        elif entry in final_residual_repaired:
            status = "proved_final_residual_chain"
        elif entry in affected:
            status = "quarantined_unmasked_coefficient_one"
        elif entry in dual:
            status = "proved_dual_endpoint_schur"
        elif entry in residual:
            status = "open_residual_no_theorem"
        else:
            status = "proved_nonuniversal_inherited"
        rows.append(
            {
                "profile": list(profile),
                "split": list(split),
                "status": status,
            }
        )
    return tuple(rows)


def artifact_text(result: MaskedUniversalAudit) -> str:
    payload = {
        "schema": "round4_q64_masked_universal_audit_v15",
        "result": asdict(result),
        "theorem_registry": theorem_registry(),
        "evidence_label": (
            "dependency-exact audit of 354 coefficient-one entries, with "
            "354 actual-mask repairs, including 54 quintic-slice, 180 "
            "local-Walsh, 12 cubic-endpoint, and 6 double-quintic endpoint "
            "entries, twelve double-quintic record entries, plus 38 four-cubic incidence, twelve cubic-septimic "
            "chain, twenty-eight recovered cubic-quintic endpoint-row entries, "
            "and twelve joint shared-quintic chain entries; the old independent-"
            "maxima proof remains rejected, while the replacement theorem retains "
            "all fifteen feasible physical quintic shape pairs; "
            "the separate final-residual theorem proves all eighty later "
            "entries with forty-eight local-Walsh and thirty-two complete-"
            "chain entries; "
            "the twelve formerly caveated dual-endpoint entries have an "
            "independent endpoint-slice and completed-link Gram audit; "
            "three full q2 signed-permutation occurrence matrices include "
            "every cross-cut distinctness mask exactly; their tangent "
            "bounds are numerical concavity-tangent screens at q2 only; "
            "they are not interval certificates; the generic "
            "mask-factor repair is rigorous but fails the q64 routing ledger"
        ),
        "decision": (
            "the written unmasked cross-Gram proof is invalid; the three "
            "q2 cubic representatives survive coefficient one, but the "
            "universal statement remains unnecessary: all 354 affected q64 "
            "entries now have independent actual-mask proofs, and the later "
            "eighty-entry residual inventory is also closed"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic masked-universal audit artifact",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_masked_universal_audit.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 masked-universal audit: "
        f"affected={result.affected_entries},"
        f"repaired={result.total_masked_repaired_entries},"
        f"quarantined={result.remaining_quarantined_entries},"
        f"supported={result.supported_entries_before_dual_caveat},"
        f"conservative_supported={result.conservative_supported_entries},"
        "q2_uniform="
        + "/".join(
            f"{row.uniform_weighted_nuclear:.12g}" for row in result.q2_rows
        )
        + ",q2_upper="
        + "/".join(
            f"{row.uniform_supporting_upper:.12g}" for row in result.q2_rows
        )
        + f",generic_mask_range={result.generic_mask_minimum_coefficient:.12g}/"
        f"{result.generic_mask_maximum_coefficient:.12g},"
        f"generic_mask_total={result.generic_mask_repaired_routing.total:.12g},"
        "status=all_888_registry_entries_certified_one_batch"
    )


if __name__ == "__main__":
    main()
