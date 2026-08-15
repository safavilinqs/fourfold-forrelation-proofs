#!/usr/bin/env python3
"""Dependency-exact, outward-rounded complete q=64 one-batch ledger."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from functools import lru_cache
from json import dumps
from math import ceil, sqrt
from pathlib import Path

from q64_dual_endpoint_schur_insertion import (
    cubic_fixed_pair_energy,
    dual_endpoint_entries,
    local_slice_coefficients,
    quintic_fixed_triple_energy,
)
from q64_final_residual_chain_contraction import (
    coefficient_map as final_residual_coefficients,
)
from q64_masked_universal_audit import (
    quarantined_coefficients,
    theorem_registry,
)
from q64_paper_target_gate import (
    DIMENSION,
    MODES,
    ORDER,
    balanced_open_entries,
    order_context,
)
import occupation_compatible_sector_optimization as occupation


ROOT = Path(__file__).resolve().parents[1]
BETA = Fraction(19, 25)
COEFFICIENT_GRID = 10**9
PROMISE_EXP_STEPS = 4096
RESERVE_THRESHOLD = Fraction(1, 3) - Fraction(1, 1000)

Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class CompleteOutwardLedger:
    order: int
    dimension: int
    sign_modes: int
    high_sector_profile_splits: int
    certified_balanced_high_sector_coefficients: int
    balanced_entries: int
    excluded_unbalanced_high_sector_entries: int
    excluded_unbalanced_high_sector_incidence_records: int
    excluded_unbalanced_high_sector_undirected_edges: int
    supported_balanced_entries: int
    open_balanced_entries: int
    dual_endpoint_certified_entries: int
    final_residual_certified_entries: int
    coefficient_grid_denominator: int
    maximum_coefficient_inflation: float
    beta_numerator: int
    beta_denominator: int
    collatz_perron_upper: str
    promise_exponent_lower: str
    promise_loss_upper: str
    total_upper: str
    reserve_threshold: str
    reserve_margin_lower: str
    passes_reserve_gate: bool
    candidate_minimum_coordinate: str
    candidate_maximum_coordinate: str
    status_counts: tuple[tuple[str, int], ...]


def strict_grid_upper(value: float) -> Fraction:
    """Place a nonnegative binary64 coefficient strictly below a rational grid point."""

    if value < 0:
        raise ValueError(("negative coefficient", value))
    if value == 0:
        return Fraction(0)
    exact_float = Fraction.from_float(value)
    units = ceil(exact_float * COEFFICIENT_GRID) + 1
    return Fraction(units, COEFFICIENT_GRID)


def complete_coefficients() -> dict[ProfileSplit, float]:
    """Assemble the certified balanced high-sector coefficients.

    The upstream routing maps also contain provisional values for unbalanced
    profile splits.  Those values were useful during theorem search but are
    not proved coefficient bounds.  They must never enter this certificate.
    For probes block diagonal in total signal photon number, only balanced
    splits occur, so the dependency-exact theorem map is precisely the 888
    entries returned here.
    """

    routed = quarantined_coefficients()
    routed.update(final_residual_coefficients())
    dual = local_slice_coefficients()[1]
    for entry in dual_endpoint_entries():
        routed[entry] = dual
    expected = set(balanced_open_entries())
    result = {entry: routed[entry] for entry in expected}
    if set(result) != expected:
        raise AssertionError(
            ("incomplete balanced high-sector coefficient map", len(result), len(expected))
        )
    if any(2 * sum(split) != sum(profile) for profile, split in result):
        raise AssertionError("unbalanced routing placeholder entered theorem map")
    return result


def coefficient_upper_map() -> dict[ProfileSplit, Fraction]:
    return {
        entry: strict_grid_upper(value)
        for entry, value in complete_coefficients().items()
    }


def _fraction_decimal_upper(value: Fraction) -> Decimal:
    with localcontext() as context:
        context.prec = 90
        context.rounding = ROUND_CEILING
        return Decimal(value.numerator) / Decimal(value.denominator)


def _fraction_decimal_lower(value: Fraction) -> Decimal:
    with localcontext() as context:
        context.prec = 90
        context.rounding = ROUND_FLOOR
        return Decimal(value.numerator) / Decimal(value.denominator)


@lru_cache(maxsize=None)
def _sqrt_integer_upper(value: int) -> Decimal:
    with localcontext() as context:
        context.prec = 90
        context.rounding = ROUND_CEILING
        root = Decimal(value).sqrt().next_plus(context)
        while root * root < value:
            root = root.next_plus(context)
        return root


def _all_profiles(mapped: dict[ProfileSplit, Fraction]) -> tuple[Profile, ...]:
    profiles = (
        (occupation.MINIMAL_PROFILE,)
        + occupation.DEGREE_SIX_PROFILES
        + occupation.DEGREE_EIGHT_PROFILES
        + occupation.KNOWN_HIGH_DEGREE_PROFILES
    )
    mapped_profiles = tuple(sorted({profile for profile, _ in mapped}))
    return profiles + tuple(profile for profile in mapped_profiles if profile not in profiles)


def all_coefficient_uppers(
    mapped: dict[ProfileSplit, Fraction],
) -> dict[ProfileSplit, Fraction]:
    """Include inherited degree-four through known-high coefficients."""

    result = dict(mapped)
    mapped_profiles = {profile for profile, _ in mapped}
    with order_context(ORDER):
        for profile in _all_profiles(mapped):
            if profile in mapped_profiles:
                continue
            for split in occupation.profile_splits(profile):
                if 2 * sum(split) != sum(profile):
                    continue
                value = occupation.coefficient(profile, split, True, None)
                result[profile, split] = strict_grid_upper(value)
    return result


def _candidate_vector(raw: dict[ProfileSplit, float]) -> tuple[Decimal, ...]:
    """Return a strictly positive floating-Perron vector used only by Collatz."""

    with order_context(ORDER):
        certificate = occupation.certificate(
            beta=float(BETA), profile_split_coefficients=raw
        )
    vector = tuple(Decimal(str(sqrt(weight))) for _, weight in certificate.occupation_weights)
    if any(value <= 0 for value in vector):
        raise AssertionError("Collatz candidate must be strictly positive")
    return vector


def collatz_perron_upper(
    all_coefficients: dict[ProfileSplit, Fraction],
    candidate: tuple[Decimal, ...],
) -> Decimal:
    """Reconstruct the 210-state matrix and apply outward Decimal Collatz."""

    states = occupation.occupation_states()
    state_index = {state: index for index, state in enumerate(states)}
    matrix = [
        [Decimal(0) for _ in states]
        for _ in states
    ]
    with localcontext() as context:
        context.prec = 80
        context.rounding = ROUND_CEILING
        for profile in _all_profiles(all_coefficients):
            attenuation = BETA ** sum(profile)
            for split in occupation.profile_splits(profile):
                if 2 * sum(split) != sum(profile):
                    continue
                coefficient = all_coefficients[profile, split]
                if not coefficient:
                    continue
                local = _fraction_decimal_upper(coefficient * attenuation)
                complement = tuple(
                    degree - selected
                    for degree, selected in zip(profile, split, strict=True)
                )
                for left_index, state in enumerate(states):
                    if any(
                        count < selected
                        for count, selected in zip(state, split, strict=True)
                    ):
                        continue
                    partner = occupation.paired_state(state, profile, split)
                    right_index = state_index.get(partner)
                    if right_index is None:
                        continue
                    left_count = occupation.multiplicity(state, split)
                    right_count = occupation.multiplicity(partner, complement)
                    if not left_count or not right_count:
                        continue
                    constant = local * _sqrt_integer_upper(left_count * right_count)
                    half = constant / Decimal(2)
                    matrix[left_index][right_index] += half
                    matrix[right_index][left_index] += half

        ratios = []
        for row, denominator in zip(matrix, candidate, strict=True):
            numerator = sum(
                (entry * value for entry, value in zip(row, candidate, strict=True)),
                Decimal(0),
            )
            ratios.append(numerator / denominator)
        return max(ratios).next_plus(context)


def promise_upper() -> tuple[Fraction, Fraction]:
    """Return a rational lower exponent and upper two-hypothesis loss.

    The biased-Rademacher proxy is at most one. Therefore the inherited
    four-chain proxy is at most ``(1+beta^2)(1+beta^4)/N``. For ``x`` equal
    to the resulting Chernoff exponent, ``exp(-x) <= (1+x/m)^(-m)``.
    """

    gap = BETA**4 - Fraction(1, 4)
    exponent = (
        gap**2
        * DIMENSION
        / (2 * (1 + BETA**2) * (1 + BETA**4))
    )
    base = 1 + exponent / PROMISE_EXP_STEPS
    two_hypothesis = Fraction(2, 1) / base**PROMISE_EXP_STEPS
    return exponent, two_hypothesis


def registry_rows() -> tuple[dict[str, object], ...]:
    raw = complete_coefficients()
    upper = coefficient_upper_map()
    statuses = {
        (tuple(row["profile"]), tuple(row["split"])): str(row["status"])
        for row in theorem_registry()
    }
    rows = []
    for entry in sorted(balanced_open_entries()):
        profile, split = entry
        status = statuses[entry]
        rows.append(
            {
                "profile": list(profile),
                "split": list(split),
                "status": status,
                "source_coefficient": raw[entry],
                "outward_upper_numerator": upper[entry].numerator,
                "outward_upper_denominator": upper[entry].denominator,
            }
        )
    return tuple(rows)


def excluded_unbalanced_high_sector_graph_counts() -> tuple[int, int]:
    """Count the formerly misclassified unbalanced high-sector incidences.

    The count is an audit witness, not a theorem input.  It demonstrates why
    the unbalanced entries cannot be called irrelevant for the earlier broad
    model: 272 profile-split/state incidences (136 undirected occupation
    edges) connect states of different total photon number within the
    at-most-six state space.
    """

    states = occupation.occupation_states()
    state_index = {state: index for index, state in enumerate(states)}
    incidences = 0
    edges: set[tuple[int, int]] = set()
    for profile in (
        profile
        for profile in occupation.HIGH_DEGREE_PROFILES
        if profile not in occupation.KNOWN_HIGH_DEGREE_PROFILES
    ):
        for split in occupation.profile_splits(profile):
            if 2 * sum(split) == sum(profile):
                continue
            for left_index, state in enumerate(states):
                if any(
                    count < selected
                    for count, selected in zip(state, split, strict=True)
                ):
                    continue
                right_index = state_index.get(
                    occupation.paired_state(state, profile, split)
                )
                if right_index is None:
                    continue
                incidences += 1
                edges.add((min(left_index, right_index), max(left_index, right_index)))
    return incidences, len(edges)


def diagnostic(
    candidate: tuple[Decimal, ...] | None = None,
) -> CompleteOutwardLedger:
    raw = complete_coefficients()
    mapped_upper = coefficient_upper_map()
    all_upper = all_coefficient_uppers(mapped_upper)
    if candidate is None:
        candidate = _candidate_vector(
            {entry: float(value) for entry, value in all_upper.items()}
        )
    if len(candidate) != len(occupation.occupation_states()):
        raise AssertionError("Collatz candidate dimension changed")
    if any(value <= 0 for value in candidate):
        raise AssertionError("Collatz candidate must be strictly positive")
    perron = collatz_perron_upper(all_upper, candidate)
    exponent, promise_fraction = promise_upper()
    promise = _fraction_decimal_upper(promise_fraction)
    with localcontext() as context:
        context.prec = 80
        context.rounding = ROUND_CEILING
        total = (perron + promise).next_plus(context)
        threshold = _fraction_decimal_lower(RESERVE_THRESHOLD)
        margin = (threshold - total).next_minus(context)
    status_counts = Counter(row["status"] for row in theorem_registry())
    balanced = set(balanced_open_entries())
    if len(balanced) != 888 or not balanced.issubset(raw):
        raise AssertionError("balanced registry is incomplete")
    if set(raw) != balanced:
        raise AssertionError("theorem coefficient map contains unbalanced placeholders")
    exact_dual_squared = cubic_fixed_pair_energy() * quintic_fixed_triple_energy()
    for entry in dual_endpoint_entries():
        if mapped_upper[entry] ** 2 < exact_dual_squared:
            raise AssertionError(("dual coefficient rounded inward", entry))
    inflations = (
        float(mapped_upper[entry] - Fraction.from_float(value))
        for entry, value in raw.items()
    )
    high_sector_profile_splits = sum(
        len(occupation.profile_splits(profile))
        for profile in occupation.HIGH_DEGREE_PROFILES
        if profile not in occupation.KNOWN_HIGH_DEGREE_PROFILES
    )
    unbalanced_incidences, unbalanced_edges = (
        excluded_unbalanced_high_sector_graph_counts()
    )
    return CompleteOutwardLedger(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        high_sector_profile_splits=high_sector_profile_splits,
        certified_balanced_high_sector_coefficients=len(raw),
        balanced_entries=len(balanced),
        excluded_unbalanced_high_sector_entries=(
            high_sector_profile_splits - len(balanced)
        ),
        excluded_unbalanced_high_sector_incidence_records=unbalanced_incidences,
        excluded_unbalanced_high_sector_undirected_edges=unbalanced_edges,
        supported_balanced_entries=len(balanced),
        open_balanced_entries=0,
        dual_endpoint_certified_entries=len(dual_endpoint_entries()),
        final_residual_certified_entries=len(final_residual_coefficients()),
        coefficient_grid_denominator=COEFFICIENT_GRID,
        maximum_coefficient_inflation=max(inflations),
        beta_numerator=BETA.numerator,
        beta_denominator=BETA.denominator,
        collatz_perron_upper=str(perron),
        promise_exponent_lower=str(_fraction_decimal_lower(exponent)),
        promise_loss_upper=str(promise),
        total_upper=str(total),
        reserve_threshold=str(threshold),
        reserve_margin_lower=str(margin),
        passes_reserve_gate=total < threshold,
        candidate_minimum_coordinate=str(min(candidate)),
        candidate_maximum_coordinate=str(max(candidate)),
        status_counts=tuple(sorted(status_counts.items())),
    )


def artifact_text(
    result: CompleteOutwardLedger,
    candidate: tuple[Decimal, ...],
) -> str:
    payload = {
        "schema": "round4_q64_complete_outward_ledger_v3",
        "result": asdict(result),
        "collatz_candidate": [str(value) for value in candidate],
        "coefficient_registry": registry_rows(),
        "rounding_contract": (
            "Every accepted binary64 theorem coefficient is replaced by a strict "
            "rational upper on a 1e-9 grid, with one additional grid unit of "
            "guard. The 210-state matrix is reconstructed with 80-digit "
            "ROUND_CEILING Decimal arithmetic and certified by the committed "
            "210-coordinate positive Collatz vector. Floating optimization is "
            "used only to discover that vector, never to verify the bound. The "
            "promise loss uses a rational "
            "Kearns-Saul-proxy relaxation and a rational exponential upper."
        ),
        "number_sector_contract": (
            "The theorem assumes every fresh batch is block diagonal in total "
            "signal photon number. Hence only splits with 2|t| = |a| occur. "
            "The 5,128 unbalanced high-sector routing entries are excluded, not "
            "certified: 272 of them would create profile-split/state incidences "
            "between different-number states (136 undirected occupation edges) "
            "in the broader vacuum-coherent model."
        ),
        "scope": (
            "Complete arbitrary-correlated-diagonal one-batch q64 ledger for "
            "fresh probes block diagonal in total signal photon number only; "
            "this artifact does not establish an adaptive lift or cover "
            "coherence between different signal-number sectors."
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    mapped_upper = coefficient_upper_map()
    all_upper = all_coefficient_uppers(mapped_upper)
    candidate = _candidate_vector(
        {entry: float(value) for entry, value in all_upper.items()}
    )
    result = diagnostic(candidate)
    text = artifact_text(result, candidate)
    if arguments.output is not None:
        arguments.output.write_text(text, encoding="utf-8")
    print(
        "q64 complete outward ledger: "
        f"supported={result.supported_balanced_entries},"
        f"open={result.open_balanced_entries},"
        f"perron={result.collatz_perron_upper},"
        f"promise={result.promise_loss_upper},"
        f"total={result.total_upper},"
        f"reserve_margin={result.reserve_margin_lower},"
        f"passes_reserve={result.passes_reserve_gate}"
    )


if __name__ == "__main__":
    main()
