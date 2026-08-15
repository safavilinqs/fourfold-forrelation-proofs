#!/usr/bin/env python3
"""Exact compatible-law diagnostic for the leading unresolved frontier.

The leading Round 4 orbit families contain block degrees five and seven, so
the inherited direct ``q=2`` physical-law diagnostic cannot represent them.
This script moves the same question to ``q=4``. It uses one shared
occupation distribution, samples one common conditional law on actual
four-block probe configurations, and inserts every leading-family matrix
entry induced by that law into a single Hermitian kernel.

All signed-permutation moments are exact rational numbers. Floating point is
used only for the final finite-dimensional eigendecompositions. The selected
``q=4`` law is also embedded in the upper-left coordinates at ``q=32`` and
reevaluated exactly. The result measures simultaneous activation and
cancellation; because the conditional law is sparse and selected
diagnostically, it is not an ``N=1024`` theorem or a certified scalar-ledger
obstruction.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass, replace
from json import dumps, loads
from math import sqrt
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent / "open_problem_forr4_passive_floor_consolidation_round_3" / "searches"
)
sys.path.insert(0, str(ROUND3_SEARCHES))

from signed_permutation_link_moment import chain_moment  # noqa: E402

from finite_size_passive_ledger import (  # noqa: E402
    final_coefficients,
    optimized_total,
)


ORDER = 4
FRONTIER_ORBITS = 16
SAMPLES_PER_STATE = 6
SEARCH_TRIALS = 24
SEED = 2026071502

Profile = tuple[int, ...]
Split = tuple[int, ...]
Entry = tuple[Profile, Split]
Configuration = tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class BasisPoint:
    state: tuple[int, ...]
    configuration: Configuration
    probability: float


@dataclass(frozen=True)
class JointDiagnostic:
    order: int
    basis_dimension: int
    support_states: int
    samples_per_state: int
    search_trials: int
    activated_orbits: int
    activated_entries: int
    nonzero_matrix_entries: int
    separate_nuclear: float
    joint_nuclear: float
    cancellation_ratio: float
    attenuated_separate_nuclear: float
    attenuated_joint_nuclear: float
    attenuated_cancellation_ratio: float
    beta: float
    selected_trial: int


def frontier() -> tuple[tuple[Entry, ...], ...]:
    path = ROOT / "artifacts" / "n1024_unresolved_orbits.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = tuple(csv.DictReader(handle))
    result: list[tuple[Entry, ...]] = []
    for row in rows[:FRONTIER_ORBITS]:
        parsed = loads(row["entries"])
        result.append(
            tuple(
                (
                    tuple(int(value) for value in profile),
                    tuple(int(value) for value in split),
                )
                for profile, split in parsed
            )
        )
    if len(result) != FRONTIER_ORBITS:
        raise AssertionError(("frontier length", len(result)))
    return tuple(result)


def occupation_law() -> tuple[float, tuple[tuple[tuple[int, ...], float], ...]]:
    coefficients, _ = final_coefficients()
    _, beta, ledger, _ = optimized_total(coefficients)
    retained = tuple((state, float(weight)) for weight, state in ledger.support)
    mass = sum(weight for _, weight in retained)
    if mass < 1 - 1e-5:
        raise AssertionError(("discarded occupation mass", 1 - mass))
    return beta, tuple((state, weight / mass) for state, weight in retained)


def sample_basis(
    rng: np.random.Generator,
    law: tuple[tuple[tuple[int, ...], float], ...],
) -> tuple[BasisPoint, ...]:
    dimension = ORDER * ORDER
    result: list[BasisPoint] = []
    for state, state_probability in law:
        seen: set[Configuration] = set()
        while len(seen) < SAMPLES_PER_STATE:
            configuration = tuple(
                tuple(
                    sorted(
                        int(value)
                        for value in rng.choice(dimension, occupation, replace=False)
                    )
                )
                for occupation in state
            )
            seen.add(configuration)
        probability = state_probability / SAMPLES_PER_STATE
        result.extend(
            BasisPoint(state, configuration, probability)
            for configuration in sorted(seen)
        )
    if abs(sum(point.probability for point in result) - 1) > 2e-14:
        raise AssertionError("basis probabilities do not sum to one")
    return tuple(result)


def differences(
    left: Configuration, right: Configuration
) -> tuple[Configuration, Profile, Split]:
    supports: list[tuple[int, ...]] = []
    split: list[int] = []
    for left_block, right_block in zip(left, right, strict=True):
        left_set = set(left_block)
        right_set = set(right_block)
        supports.append(tuple(sorted(left_set ^ right_set)))
        split.append(len(left_set - right_set))
    profile = tuple(len(support) for support in supports)
    return tuple(supports), profile, tuple(split)


def matrices(
    basis: tuple[BasisPoint, ...],
    orbit_entries: tuple[tuple[Entry, ...], ...],
    order: int,
) -> tuple[tuple[np.ndarray, ...], int, int, int]:
    entry_to_orbit = {
        entry: orbit for orbit, entries in enumerate(orbit_entries) for entry in entries
    }
    dimension = len(basis)
    result = tuple(np.zeros((dimension, dimension), dtype=float) for _ in orbit_entries)
    activated_entries: set[Entry] = set()
    nonzero = 0
    for row, left in enumerate(basis):
        for column in range(row + 1, dimension):
            right = basis[column]
            supports, profile, split = differences(
                left.configuration, right.configuration
            )
            orbit = entry_to_orbit.get((profile, split))
            if orbit is None:
                continue
            exact = chain_moment(order, supports)
            if not exact:
                continue
            weighted = sqrt(left.probability * right.probability) * float(exact)
            result[orbit][row, column] = weighted
            result[orbit][column, row] = weighted
            activated_entries.add((profile, split))
            nonzero += 2
    activated_orbits = sum(np.any(matrix) for matrix in result)
    return result, int(activated_orbits), len(activated_entries), nonzero


def nuclear(matrix: np.ndarray) -> float:
    return float(np.abs(np.linalg.eigvalsh(matrix)).sum())


def score(
    basis: tuple[BasisPoint, ...],
    orbit_entries: tuple[tuple[Entry, ...], ...],
    beta: float,
    order: int = ORDER,
) -> tuple[tuple[int, int, float], JointDiagnostic]:
    orbit_matrices, activated, entries, nonzero = matrices(basis, orbit_entries, order)
    separate_values = tuple(nuclear(matrix) for matrix in orbit_matrices)
    separate = sum(separate_values)
    joint = nuclear(sum(orbit_matrices, start=np.zeros_like(orbit_matrices[0])))
    attenuated_matrices = tuple(
        matrix * beta ** sum(orbit_entries[index][0][0])
        for index, matrix in enumerate(orbit_matrices)
    )
    attenuated_values = tuple(nuclear(matrix) for matrix in attenuated_matrices)
    attenuated_separate = sum(attenuated_values)
    attenuated_joint = nuclear(
        sum(attenuated_matrices, start=np.zeros_like(attenuated_matrices[0]))
    )
    result = JointDiagnostic(
        order=order,
        basis_dimension=len(basis),
        support_states=len({point.state for point in basis}),
        samples_per_state=SAMPLES_PER_STATE,
        search_trials=SEARCH_TRIALS,
        activated_orbits=activated,
        activated_entries=entries,
        nonzero_matrix_entries=nonzero,
        separate_nuclear=separate,
        joint_nuclear=joint,
        cancellation_ratio=joint / separate if separate else 0.0,
        attenuated_separate_nuclear=attenuated_separate,
        attenuated_joint_nuclear=attenuated_joint,
        attenuated_cancellation_ratio=(
            attenuated_joint / attenuated_separate if attenuated_separate else 0.0
        ),
        beta=beta,
        selected_trial=-1,
    )
    # First maximise simultaneous family activation, then the number of
    # realised split entries, then the actual attenuated joint norm.
    return (activated, entries, attenuated_joint), result


def embed_configuration(
    configuration: Configuration, target_order: int
) -> Configuration:
    if target_order < ORDER:
        raise ValueError(("embedding order", ORDER, target_order))
    return tuple(
        tuple(
            (coordinate // ORDER) * target_order + coordinate % ORDER
            for coordinate in block
        )
        for block in configuration
    )


def embed_basis(
    basis: tuple[BasisPoint, ...], target_order: int
) -> tuple[BasisPoint, ...]:
    return tuple(
        BasisPoint(
            point.state,
            embed_configuration(point.configuration, target_order),
            point.probability,
        )
        for point in basis
    )


def diagnostics() -> tuple[JointDiagnostic, JointDiagnostic]:
    orbit_entries = frontier()
    beta, law = occupation_law()
    rng = np.random.default_rng(SEED)
    best_key: tuple[int, int, float] | None = None
    best: JointDiagnostic | None = None
    best_basis: tuple[BasisPoint, ...] | None = None
    for trial in range(SEARCH_TRIALS):
        basis = sample_basis(rng, law)
        key, candidate = score(basis, orbit_entries, beta, ORDER)
        if best_key is None or key > best_key:
            best_key = key
            best = replace(candidate, selected_trial=trial)
            best_basis = basis
    if best is None or best_basis is None:
        raise AssertionError("joint-law search produced no candidate")
    _, embedded = score(embed_basis(best_basis, 32), orbit_entries, beta, 32)
    return best, replace(embedded, selected_trial=best.selected_trial)


def diagnostic() -> JointDiagnostic:
    """Return the searched q=4 diagnostic for compatibility with callers."""

    return diagnostics()[0]


def format_result(result: JointDiagnostic) -> str:
    label = "q4" if result.order == 4 else "embedded_q32"
    return (
        f"{label}:basis={result.basis_dimension},"
        f"states={result.support_states},"
        f"samples_per_state={result.samples_per_state},"
        f"trials={result.search_trials},"
        f"selected_trial={result.selected_trial},"
        f"activated_orbits={result.activated_orbits}/{FRONTIER_ORBITS},"
        f"activated_entries={result.activated_entries},"
        f"nonzero_entries={result.nonzero_matrix_entries},"
        f"separate={result.separate_nuclear:.12g},"
        f"joint={result.joint_nuclear:.12g},"
        f"ratio={result.cancellation_ratio:.12g},"
        f"attenuated_separate={result.attenuated_separate_nuclear:.12g},"
        f"attenuated_joint={result.attenuated_joint_nuclear:.12g},"
        f"attenuated_ratio={result.attenuated_cancellation_ratio:.12g}"
    )


def artifact_text(q4: JointDiagnostic, q32: JointDiagnostic) -> str:
    payload = {
        "schema": "round4_joint_impact_sparse_diagnostic_v1",
        "frontier_orbits": FRONTIER_ORBITS,
        "selection_rule": (
            "maximize activated orbits, then activated entries, then q4 "
            "attenuated joint nuclear norm"
        ),
        "q4": asdict(q4),
        "embedded_q32": asdict(q32),
        "evidence_label": (
            "exact plant moments with floating eigendecomposition; "
            "compatible physical-law diagnostic, not theorem coefficient"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic JSON result under artifacts/",
    )
    arguments = parser.parse_args()
    q4, q32 = diagnostics()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "joint_impact_sparse_diagnostic.json"
        path.write_text(artifact_text(q4, q32), encoding="utf-8")
    print(
        "joint q4 impact-frontier diagnostic: "
        f"{format_result(q4)};{format_result(q32)};"
        f"beta={q4.beta:.12g},"
        "status=diagnostic_not_N1024_certificate"
    )


if __name__ == "__main__":
    main()
