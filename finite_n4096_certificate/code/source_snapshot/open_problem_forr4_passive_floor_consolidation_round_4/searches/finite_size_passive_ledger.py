#!/usr/bin/env python3
"""Generate the authoritative Round 4 finite-size routing ledger.

This audit does not turn the current ``1/q`` completion into a theorem.  It
reconstructs the calculation from the Round 3 source modules, labels every
balanced high-degree entry by the evidence actually available, and records
the exact point at which the current signed-permutation route remains open.

The signed-permutation implementation uses ``N=q^2`` with power-of-two
``q``.  Consequently only ``N=256,1024,4096`` in the declared Round 4 window
have the current witness geometry.  The complete coefficient ledger is
currently calibrated only at ``N=1024``.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from json import dumps, loads
from math import isqrt
from pathlib import Path
import sys
from typing import Callable

from scipy.optimize import brentq, minimize_scalar


ROOT = Path(__file__).resolve().parents[1]
ROUND3 = ROOT.parent / "open_problem_forr4_passive_floor_consolidation_round_3"
ROUND3_SEARCHES = ROUND3 / "searches"
sys.path.insert(0, str(ROUND3_SEARCHES))

from adjacent_balanced_cubic_slice_contraction import (  # noqa: E402
    adjacent_balanced_coefficient,
    adjacent_balanced_orbit_entries,
)
from adjacent_balanced_row_slice_contraction import (  # noqa: E402
    adjacent_balanced_row_coefficient,
    target_orbit_entries as adjacent_balanced_row_orbit_entries,
)
from attenuation_promise_concentration import (  # noqa: E402
    extended_euclidean_promise_concentration,
    promise_concentration,
)
from column_cubic_quintic_row_contraction import (  # noqa: E402
    column_cubic_quintic_coefficient,
    column_cubic_quintic_orbit_entries,
)
from double_endpoint_cubic_quintic_row_contraction import (  # noqa: E402
    double_endpoint_cubic_quintic_coefficient,
    double_endpoint_cubic_quintic_orbit_entries,
)
from high_degree_record_incidence_frontier import (  # noqa: E402
    certified_coefficients,
    compatible_record_triples,
    dose_six_relevant_entries,
    split_perron_sensitivity,
    symmetry_orbits,
)
from internal_singleton_shared_law_contraction import (  # noqa: E402
    internal_singleton_coefficient,
    internal_singleton_orbit_entries,
)
from leading_balanced_disjointness_contraction import (  # noqa: E402
    leading_balanced_coefficient,
    leading_balanced_orbit_entries,
)
from middle_cubic_quintic_pair_contraction import (  # noqa: E402
    middle_cubic_quintic_pair_coefficient,
    middle_cubic_quintic_pair_orbit_entries,
)
from occupation_compatible_sector_optimization import (  # noqa: E402
    BETA,
    certificate,
)
from repaired_open_profile_budget import (  # noqa: E402
    ADJACENT_SPLIT_ENTRIES,
    LEADING_SPLIT_ENTRIES,
    coarse_open_completion_coefficients,
    coarse_open_completion_target,
    forced_coefficients,
)
from separated_balanced_endpoint_slice_contraction import (  # noqa: E402
    separated_balanced_coefficient,
    separated_balanced_orbit_entries,
)
from whole_cubic_middle_pair_contraction import (  # noqa: E402
    whole_cubic_middle_pair_coefficient,
    whole_cubic_middle_pair_orbit_entries,
)
from whole_cubic_quintic_triple_contraction import (  # noqa: E402
    whole_cubic_quintic_triple_coefficient,
    whole_cubic_quintic_triple_orbit_entries,
)


DIMENSION = 1024
ORDER = 32
THRESHOLD = 1 / 3
DECLARED_WINDOW = (256, 512, 1024, 2048, 4096)
ROUND3_NOTES = (
    "../open_problem_forr4_passive_floor_consolidation_round_3/notes"
)

Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


def round3_note(filename: str) -> str:
    """Return a Round-4-root-relative path to inherited evidence."""

    return f"{ROUND3_NOTES}/{filename}"


@dataclass(frozen=True)
class TheoremFamily:
    """One proved complement/reversal orbit used by the final diagnostic."""

    name: str
    entries: tuple[ProfileSplit, ...]
    coefficient: float
    source: str


@dataclass(frozen=True)
class LedgerTotals:
    """Reconstructed numerical checkpoints and their evidence labels."""

    degree_eight_total: float
    known_high_total: float
    coarse_before_ten_total: float
    ten_theorem_beta: float
    ten_theorem_perron_upper: float
    ten_theorem_promise_loss: float
    ten_theorem_total: float
    ten_theorem_margin: float


@dataclass(frozen=True)
class LedgerAudit:
    """Complete in-memory audit returned to tests and artifact writers."""

    totals: LedgerTotals
    rows: tuple[dict[str, object], ...]
    unresolved_orbits: tuple[dict[str, object], ...]
    window: tuple[dict[str, object], ...]
    theorem_families: tuple[TheoremFamily, ...]
    counts: dict[str, int]
    impact_frontier: dict[str, float | int]
    top_unresolved_gate: float
    top_provisional_gate: float


def theorem_families(dimension: int = DIMENSION) -> tuple[TheoremFamily, ...]:
    """Return the ten accepted arbitrary-law orbit theorems."""

    specifications: tuple[
        tuple[
            str,
            Callable[[], tuple[ProfileSplit, ...]],
            Callable[[int], float],
            str,
        ],
        ...,
    ] = (
        (
            "leading_disjointness",
            leading_balanced_orbit_entries,
            leading_balanced_coefficient,
            round3_note("leading_balanced_disjointness_contraction.md"),
        ),
        (
            "adjacent_cubic_slice",
            adjacent_balanced_orbit_entries,
            adjacent_balanced_coefficient,
            round3_note("adjacent_balanced_cubic_slice_contraction.md"),
        ),
        (
            "separated_endpoint_slice",
            separated_balanced_orbit_entries,
            separated_balanced_coefficient,
            round3_note("separated_balanced_endpoint_slice_contraction.md"),
        ),
        (
            "internal_singleton_shared_law",
            internal_singleton_orbit_entries,
            internal_singleton_coefficient,
            round3_note("internal_singleton_shared_law_contraction.md"),
        ),
        (
            "column_cubic_quintic_row",
            column_cubic_quintic_orbit_entries,
            column_cubic_quintic_coefficient,
            round3_note("column_cubic_quintic_row_contraction.md"),
        ),
        (
            "adjacent_balanced_row",
            adjacent_balanced_row_orbit_entries,
            adjacent_balanced_row_coefficient,
            round3_note("adjacent_balanced_row_slice_contraction.md"),
        ),
        (
            "whole_cubic_quintic_triple",
            whole_cubic_quintic_triple_orbit_entries,
            whole_cubic_quintic_triple_coefficient,
            round3_note("whole_cubic_quintic_triple_contraction.md"),
        ),
        (
            "middle_cubic_quintic_pair",
            middle_cubic_quintic_pair_orbit_entries,
            middle_cubic_quintic_pair_coefficient,
            round3_note("middle_cubic_quintic_pair_contraction.md"),
        ),
        (
            "whole_cubic_middle_pair",
            whole_cubic_middle_pair_orbit_entries,
            whole_cubic_middle_pair_coefficient,
            round3_note("whole_cubic_middle_pair_contraction.md"),
        ),
        (
            "double_endpoint_cubic_quintic_row",
            double_endpoint_cubic_quintic_orbit_entries,
            double_endpoint_cubic_quintic_coefficient,
            round3_note("double_endpoint_cubic_quintic_row_contraction.md"),
        ),
    )
    return tuple(
        TheoremFamily(
            name=name,
            entries=entries(),
            coefficient=coefficient(dimension),
            source=source,
        )
        for name, entries, coefficient, source in specifications
    )


def final_coefficients(
    dimension: int = DIMENSION,
) -> tuple[dict[ProfileSplit, float], tuple[TheoremFamily, ...]]:
    """Apply the ten theorems to the inherited diagnostic coefficient map."""

    if dimension != DIMENSION:
        raise ValueError(("complete ledger calibrated only at N=1024", dimension))
    coefficients = coarse_open_completion_coefficients()
    families = theorem_families(dimension)
    seen: set[ProfileSplit] = set()
    for family in families:
        if len(family.entries) != 4 or seen.intersection(family.entries):
            raise AssertionError(("theorem orbit overlap", family))
        seen.update(family.entries)
        for entry in family.entries:
            coefficients[entry] = family.coefficient
    return coefficients, families


def optimized_total(
    coefficients: dict[ProfileSplit, float],
    dimension: int = DIMENSION,
) -> tuple[float, float, object, object]:
    """Optimize attenuation and return the diagnostic ledger components."""

    def total(beta: float) -> float:
        ledger = certificate(
            beta=beta,
            profile_split_coefficients=coefficients,
        )
        promise = extended_euclidean_promise_concentration(dimension, beta)
        return ledger.supporting_upper + promise.two_hypothesis_loss

    optimum = minimize_scalar(
        total,
        bounds=(0.75, 0.81),
        method="bounded",
        options={"xatol": 1e-13},
    )
    beta = float(optimum.x)
    ledger = certificate(
        beta=beta,
        profile_split_coefficients=coefficients,
    )
    promise = extended_euclidean_promise_concentration(dimension, beta)
    return float(optimum.fun), beta, ledger, promise


def _status_maps(
    families: tuple[TheoremFamily, ...],
) -> tuple[dict[ProfileSplit, TheoremFamily], set[ProfileSplit], set[ProfileSplit]]:
    theorem_by_entry = {
        entry: family for family in families for entry in family.entries
    }
    physical_lower = set(forced_coefficients()) | set(LEADING_SPLIT_ENTRIES)
    local_slice = set(ADJACENT_SPLIT_ENTRIES)
    if physical_lower.intersection(local_slice):
        raise AssertionError("diagnostic substitution orbits overlap")
    return theorem_by_entry, physical_lower, local_slice


def _entry_status(
    entry: ProfileSplit,
    theorem_by_entry: dict[ProfileSplit, TheoremFamily],
    physical_lower: set[ProfileSplit],
    local_slice: set[ProfileSplit],
) -> tuple[str, str, str]:
    if entry in theorem_by_entry:
        family = theorem_by_entry[entry]
        return "proved_arbitrary_upper", family.name, family.source
    if entry in physical_lower:
        return (
            "physical_lower_witness_substitution",
            "diagnostic_physical_orbit",
            round3_note("repaired_open_profile_budget.md"),
        )
    if entry in local_slice:
        return (
            "local_slice_substitution_not_arbitrary_upper",
            "diagnostic_adjacent_slice",
            round3_note("adjacent_cubic_quintic_record_gate.md"),
        )
    return (
        "provisional_common_1_over_q_target",
        "unresolved_balanced_frontier",
        round3_note("high_degree_record_incidence_frontier.md"),
    )


def _tuple_json(value: tuple[int, ...]) -> str:
    return dumps(value, separators=(",", ":"))


def _profile_family(profile: Profile) -> Profile:
    reverse = tuple(reversed(profile))
    return min(profile, reverse)


def _window_rows(
    theorem_specs: tuple[TheoremFamily, ...],
    ten_total: float,
    ten_margin: float,
) -> tuple[dict[str, object], ...]:
    result = []
    for dimension in DECLARED_WINDOW:
        order = isqrt(dimension)
        if order * order != dimension:
            status = "unsupported_by_current_N_equals_q_squared_witness"
            reason = "sqrt(N) is not an integer"
            order_value: int | None = None
        elif order & (order - 1):
            status = "unsupported_non_power_of_two_order"
            reason = "the signed-permutation Walsh construction requires power-of-two q"
            order_value = order
        elif dimension == DIMENSION:
            status = "diagnostic_only_complete_q32_ledger"
            reason = "848 balanced entries still lack arbitrary-law upper bounds"
            order_value = order
        else:
            status = "witness_geometry_supported_but_full_ledger_not_calibrated"
            reason = "accepted sectors and physical diagnostic orbits remain q=32-specific"
            order_value = order
        row: dict[str, object] = {
            "N": dimension,
            "q": order_value,
            "M": 4 * dimension,
            "status": status,
            "reason": reason,
            "ten_theorem_total": "",
            "margin_to_one_third": "",
        }
        if dimension == DIMENSION:
            row["ten_theorem_total"] = ten_total
            row["margin_to_one_third"] = ten_margin
        elif order_value is not None:
            coefficients = [
                family.coefficient
                for family in theorem_families(dimension)
            ]
            row["proved_orbit_coefficient_min"] = min(coefficients)
            row["proved_orbit_coefficient_max"] = max(coefficients)
        result.append(row)
    if tuple(family.name for family in theorem_specs) != tuple(
        family.name for family in theorem_families()
    ):
        raise AssertionError("theorem family order changed")
    return tuple(result)


def build_audit() -> LedgerAudit:
    """Reconstruct every inherited checkpoint and classify all 888 entries."""

    coefficients, families = final_coefficients()
    total, beta, ledger, promise = optimized_total(coefficients)

    degree_eight_beta = 313 / 400
    degree_eight = certificate(
        beta=degree_eight_beta,
        include_known_high_degree=False,
    )
    degree_eight_total = (
        degree_eight.supporting_upper
        + promise_concentration(DIMENSION, degree_eight_beta).two_hypothesis_loss
    )
    known = certificate(beta=BETA)
    known_total = (
        known.supporting_upper
        + promise_concentration(DIMENSION, BETA).two_hypothesis_loss
    )
    coarse_total = coarse_open_completion_target().optimized_total
    totals = LedgerTotals(
        degree_eight_total=degree_eight_total,
        known_high_total=known_total,
        coarse_before_ten_total=coarse_total,
        ten_theorem_beta=beta,
        ten_theorem_perron_upper=ledger.supporting_upper,
        ten_theorem_promise_loss=promise.two_hypothesis_loss,
        ten_theorem_total=total,
        ten_theorem_margin=THRESHOLD - total,
    )

    relevant = dose_six_relevant_entries()
    relevant_set = set(relevant)
    theorem_by_entry, physical_lower, local_slice = _status_maps(families)
    if not set(theorem_by_entry) <= relevant_set:
        raise AssertionError("a ten-theorem entry left the balanced shell")
    if not (physical_lower | local_slice) <= relevant_set:
        raise AssertionError("a diagnostic substitution left the balanced shell")

    weights = dict(ledger.occupation_weights)
    rows = []
    for profile, split in relevant:
        entry = (profile, split)
        status, family, source = _entry_status(
            entry,
            theorem_by_entry,
            physical_lower,
            local_slice,
        )
        sensitivity = split_perron_sensitivity(
            profile,
            split,
            beta,
            weights,
        )
        coefficient = coefficients[entry]
        rows.append(
            {
                "profile": _tuple_json(profile),
                "split": _tuple_json(split),
                "degree": sum(profile),
                "selected_degree": sum(split),
                "profile_reversal_family": _tuple_json(_profile_family(profile)),
                "compatible_record_triples": dumps(
                    compatible_record_triples(profile),
                    separators=(",", ":"),
                ),
                "coefficient": coefficient,
                "coefficient_times_q": coefficient * ORDER,
                "evidence_status": status,
                "family": family,
                "source": source,
                "perron_sensitivity": sensitivity,
                "perron_contribution": sensitivity * coefficient,
            }
        )

    unresolved_statuses = {
        "physical_lower_witness_substitution",
        "local_slice_substitution_not_arbitrary_upper",
        "provisional_common_1_over_q_target",
    }
    row_by_entry = {
        (tuple(loads(str(row["profile"]))), tuple(loads(str(row["split"])))): row
        for row in rows
    }
    unresolved_orbits = []
    for orbit in symmetry_orbits(relevant):
        orbit_rows = tuple(row_by_entry[entry] for entry in orbit)
        if not any(row["evidence_status"] in unresolved_statuses for row in orbit_rows):
            continue
        statuses = sorted({str(row["evidence_status"]) for row in orbit_rows})
        unresolved_orbits.append(
            {
                "entries": dumps(orbit, separators=(",", ":")),
                "size": len(orbit),
                "statuses": ";".join(statuses),
                "profile_reversal_family": orbit_rows[0]["profile_reversal_family"],
                "coefficient_min": min(float(row["coefficient"]) for row in orbit_rows),
                "coefficient_max": max(float(row["coefficient"]) for row in orbit_rows),
                "perron_sensitivity": sum(
                    float(row["perron_sensitivity"]) for row in orbit_rows
                ),
                "perron_contribution": sum(
                    float(row["perron_contribution"]) for row in orbit_rows
                ),
            }
        )
    unresolved_orbits.sort(
        key=lambda row: float(row["perron_contribution"]),
        reverse=True,
    )

    def orbit_entries(orbit_row: dict[str, object]) -> tuple[ProfileSplit, ...]:
        return tuple(
            (tuple(profile), tuple(split))
            for profile, split in loads(str(orbit_row["entries"]))
        )

    def orbit_gate(orbit_row: dict[str, object]) -> float:
        entries = orbit_entries(orbit_row)

        def total_with_value(value: float) -> float:
            trial = dict(coefficients)
            for entry in entries:
                trial[entry] = value
            return optimized_total(trial)[0]

        current = max(coefficients[entry] for entry in entries)
        return float(
            brentq(
                lambda value: total_with_value(value) - THRESHOLD,
                current,
                0.3,
                xtol=1e-14,
            )
        )

    top_gate = orbit_gate(unresolved_orbits[0])
    provisional_orbits = tuple(
        row
        for row in unresolved_orbits
        if row["statuses"] == "provisional_common_1_over_q_target"
    )
    top_provisional_gate = orbit_gate(provisional_orbits[0])

    counts = {
        "open_profile_splits": len(coarse_open_completion_coefficients()),
        "off_face_generic_certified_at_1_over_q": len(certified_coefficients()),
        "balanced_entries": len(rows),
        "proved_arbitrary_upper_entries": sum(
            row["evidence_status"] == "proved_arbitrary_upper" for row in rows
        ),
        "physical_lower_witness_substitutions": sum(
            row["evidence_status"] == "physical_lower_witness_substitution"
            for row in rows
        ),
        "local_slice_substitutions": sum(
            row["evidence_status"]
            == "local_slice_substitution_not_arbitrary_upper"
            for row in rows
        ),
        "provisional_common_target_entries": sum(
            row["evidence_status"] == "provisional_common_1_over_q_target"
            for row in rows
        ),
        "unresolved_entries": sum(
            row["evidence_status"] in unresolved_statuses for row in rows
        ),
        "unresolved_orbits": len(unresolved_orbits),
    }
    contributions = tuple(
        float(row["perron_contribution"]) for row in unresolved_orbits
    )
    total_unresolved_contribution = sum(contributions)

    def orbits_reaching(fraction: float) -> int:
        cumulative = 0.0
        for index, contribution in enumerate(contributions, start=1):
            cumulative += contribution
            if cumulative >= fraction * total_unresolved_contribution:
                return index
        raise AssertionError(("impact frontier", fraction, cumulative))

    impact_frontier: dict[str, float | int] = {
        "total_unresolved_perron_contribution": total_unresolved_contribution,
        "leading_orbit_fraction": contributions[0]
        / total_unresolved_contribution,
        "leading_16_orbits_fraction": sum(contributions[:16])
        / total_unresolved_contribution,
        "leading_51_orbits_fraction": sum(contributions[:51])
        / total_unresolved_contribution,
        "orbits_reaching_50_percent": orbits_reaching(0.5),
        "orbits_reaching_90_percent": orbits_reaching(0.9),
    }
    window = _window_rows(families, total, THRESHOLD - total)
    return LedgerAudit(
        totals=totals,
        rows=tuple(rows),
        unresolved_orbits=tuple(unresolved_orbits),
        window=window,
        theorem_families=families,
        counts=counts,
        impact_frontier=impact_frontier,
        top_unresolved_gate=float(top_gate),
        top_provisional_gate=top_provisional_gate,
    )


def _write_csv(path: Path, rows: tuple[dict[str, object], ...]) -> None:
    if not rows:
        raise ValueError(("empty CSV artifact", path))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=tuple(rows[0]),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_artifacts(audit: LedgerAudit, output_dir: Path) -> None:
    """Write deterministic human-inspectable ledger artifacts."""

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "n1024_balanced_ledger.csv", audit.rows)
    _write_csv(output_dir / "n1024_unresolved_orbits.csv", audit.unresolved_orbits)
    _write_csv(output_dir / "finite_size_window.csv", audit.window)
    summary = {
        "label": "diagnostic_routing_ledger_not_passive_lower_bound",
        "N": DIMENSION,
        "q": ORDER,
        "M": 4 * DIMENSION,
        "counts": audit.counts,
        "impact_frontier": audit.impact_frontier,
        "totals": asdict(audit.totals),
        "theorem_families": [asdict(family) for family in audit.theorem_families],
        "top_unresolved_orbit": audit.unresolved_orbits[0],
        "top_unresolved_gate": audit.top_unresolved_gate,
        "top_provisional_orbit": next(
            row
            for row in audit.unresolved_orbits
            if row["statuses"] == "provisional_common_1_over_q_target"
        ),
        "top_provisional_gate": audit.top_provisional_gate,
        "completion_missing": [
            "848 arbitrary-law balanced-entry upper bounds",
            "outward-rounded interval Perron certificate",
            "adaptive posterior-selection lift",
            "experimental feasibility assessment",
        ],
    }
    (output_dir / "finite_size_ledger_summary.json").write_text(
        dumps(summary, indent=2, sort_keys=True) + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="write CSV and JSON artifacts to this directory",
    )
    arguments = parser.parse_args()
    audit = build_audit()
    if arguments.output_dir is not None:
        write_artifacts(audit, arguments.output_dir)
    top = audit.unresolved_orbits[0]
    print(
        "finite-size passive ledger: "
        f"balanced={audit.counts['balanced_entries']},"
        f"proved={audit.counts['proved_arbitrary_upper_entries']},"
        f"unresolved={audit.counts['unresolved_entries']},"
        f"unresolved_orbits={audit.counts['unresolved_orbits']},"
        f"beta={audit.totals.ten_theorem_beta:.15g},"
        f"total={audit.totals.ten_theorem_total:.15g},"
        f"margin={audit.totals.ten_theorem_margin:.15g},"
        f"top={top['entries']},"
        f"top_gate={audit.top_unresolved_gate:.15g},"
        "label=diagnostic_not_theorem"
    )


if __name__ == "__main__":
    main()
