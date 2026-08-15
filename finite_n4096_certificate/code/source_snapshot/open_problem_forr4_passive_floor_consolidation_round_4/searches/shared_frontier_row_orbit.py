#!/usr/bin/env python3
"""Native-q32 row-translation physical law for the shared frontier.

For every occupation state in the current Perron support, choose one base
four-block configuration and average it over the common row-XOR translations
of all four blocks.  At q=32 this is a legal 960-configuration diagonal probe
law.  The script searches deterministic base configurations for simultaneous
frontier activation, evaluates all matching signed-permutation moments
exactly, and computes separate and joint trace norms in floating point.

This is a compatible physical lower-witness diagnostic.  It is neither an
arbitrary-law upper contraction nor an interval-certified kill certificate.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from json import dumps
from math import sqrt
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent / "open_problem_forr4_passive_floor_consolidation_round_3" / "searches"
)
sys.path.insert(0, str(ROUND3_SEARCHES))

from double_endpoint_occupation_optimization import (  # noqa: E402
    occupation_states,
)

from joint_impact_sparse_q4 import occupation_law  # noqa: E402
from shared_frontier_structure import (  # noqa: E402
    compatible_edges,
    frontier_rows,
    parsed_entries,
)
from signed_permutation_link_moment import chain_moment  # noqa: E402


ORDER = 32
TRIALS = 9
SEED = 2026071600

Profile = tuple[int, ...]
Split = tuple[int, ...]
State = tuple[int, ...]
Entry = tuple[Profile, Split]
Configuration = tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class Activation:
    trial: int
    activated_orbits: int
    activated_entries: int
    matching_pairs: int
    nonzero_orbits: int
    nonzero_entries: int
    nonzero_pairs: int
    absolute_moment_sum: float


@dataclass(frozen=True)
class RowOrbitDiagnostic:
    order: int
    basis_dimension: int
    occupation_states: int
    translation_orbit_size: int
    searched_trials: int
    selected_trial: int
    relevant_occupation_edges: int
    activated_orbits: int
    activated_entries: int
    matching_pairs: int
    nonzero_pairs: int
    attenuated_separate_nuclear: float
    attenuated_joint_nuclear: float
    attenuated_cancellation_ratio: float
    beta: float
    current_frontier_perron_contribution: float
    separate_to_current_frontier_ratio: float


def frontier() -> tuple[tuple[Entry, ...], ...]:
    return tuple(parsed_entries(row) for row in frontier_rows())


def translate(configuration: Configuration, shift: int) -> Configuration:
    return tuple(
        tuple(
            sorted(
                ((coordinate // ORDER) ^ shift) * ORDER + coordinate % ORDER
                for coordinate in block
            )
        )
        for block in configuration
    )


def difference(
    left: Configuration, right: Configuration
) -> tuple[Configuration, Entry]:
    supports: list[tuple[int, ...]] = []
    profile: list[int] = []
    split: list[int] = []
    for left_block, right_block in zip(left, right, strict=True):
        left_set = set(left_block)
        right_set = set(right_block)
        support = tuple(sorted(left_set ^ right_set))
        supports.append(support)
        profile.append(len(support))
        split.append(len(left_set - right_set))
    return tuple(supports), (tuple(profile), tuple(split))


def random_bases(trial: int, states: tuple[State, ...]) -> dict[State, Configuration]:
    rng = np.random.default_rng(SEED + trial)
    dimension = ORDER * ORDER
    result: dict[State, Configuration] = {}
    for state in states:
        while True:
            candidate = tuple(
                tuple(
                    sorted(
                        int(value)
                        for value in rng.choice(dimension, occupation, replace=False)
                    )
                )
                for occupation in state
            )
            orbit = {translate(candidate, shift) for shift in range(ORDER)}
            if len(orbit) == ORDER:
                result[state] = candidate
                break
    return result


def relevant_edges(
    orbit_entries: tuple[tuple[Entry, ...], ...], states: tuple[State, ...]
) -> tuple[tuple[State, State], ...]:
    available = set(states)
    all_states = tuple(occupation_states())
    result: set[tuple[State, State]] = set()
    for entries in orbit_entries:
        for entry in entries:
            for edge in compatible_edges(entry, all_states):
                if edge[0] in available and edge[1] in available:
                    result.add(edge)
    return tuple(sorted(result))


def translated_bases(
    bases: dict[State, Configuration], states: tuple[State, ...]
) -> dict[State, tuple[Configuration, ...]]:
    return {
        state: tuple(translate(bases[state], shift) for shift in range(ORDER))
        for state in states
    }


def activation(
    trial: int,
    orbit_entries: tuple[tuple[Entry, ...], ...],
    states: tuple[State, ...],
    edges: tuple[tuple[State, State], ...],
) -> tuple[Activation, dict[State, Configuration]]:
    entry_to_orbit = {
        entry: orbit_index
        for orbit_index, entries in enumerate(orbit_entries)
        for entry in entries
    }
    bases = random_bases(trial, states)
    translated = translated_bases(bases, states)
    activated_orbits: set[int] = set()
    activated_entries: set[Entry] = set()
    nonzero_orbits: set[int] = set()
    nonzero_entries: set[Entry] = set()
    matching_pairs = 0
    nonzero_pairs = 0
    absolute_moment_sum = 0.0
    for left_state, right_state in edges:
        for left in translated[left_state]:
            for right in translated[right_state]:
                supports, entry = difference(left, right)
                orbit_index = entry_to_orbit.get(entry)
                if orbit_index is None:
                    continue
                activated_orbits.add(orbit_index)
                activated_entries.add(entry)
                matching_pairs += 1
                exact = chain_moment(ORDER, supports)
                if exact:
                    nonzero_orbits.add(orbit_index)
                    nonzero_entries.add(entry)
                    nonzero_pairs += 1
                    absolute_moment_sum += abs(float(exact))
    return (
        Activation(
            trial=trial,
            activated_orbits=len(activated_orbits),
            activated_entries=len(activated_entries),
            matching_pairs=matching_pairs,
            nonzero_orbits=len(nonzero_orbits),
            nonzero_entries=len(nonzero_entries),
            nonzero_pairs=nonzero_pairs,
            absolute_moment_sum=absolute_moment_sum,
        ),
        bases,
    )


def selected_bases(
    orbit_entries: tuple[tuple[Entry, ...], ...],
    states: tuple[State, ...],
    edges: tuple[tuple[State, State], ...],
) -> tuple[Activation, dict[State, Configuration]]:
    best: Activation | None = None
    best_bases: dict[State, Configuration] | None = None
    for trial in range(TRIALS):
        candidate, bases = activation(trial, orbit_entries, states, edges)
        key = (
            candidate.nonzero_orbits,
            candidate.nonzero_entries,
            candidate.absolute_moment_sum,
            candidate.nonzero_pairs,
        )
        best_key = (
            (
                best.nonzero_orbits,
                best.nonzero_entries,
                best.absolute_moment_sum,
                best.nonzero_pairs,
            )
            if best is not None
            else None
        )
        if best_key is None or key > best_key:
            best = candidate
            best_bases = bases
    if best is None or best_bases is None:
        raise AssertionError("row-orbit search returned no candidate")
    return best, best_bases


def nuclear_from_entries(
    entries: list[tuple[int, int, float]],
) -> float:
    active = sorted({index for left, right, _ in entries for index in (left, right)})
    local_index = {index: position for position, index in enumerate(active)}
    matrix = np.zeros((len(active), len(active)))
    for left, right, value in entries:
        matrix[local_index[left], local_index[right]] += value
        matrix[local_index[right], local_index[left]] += value
    return float(np.abs(np.linalg.eigvalsh(matrix)).sum())


def diagnostic() -> RowOrbitDiagnostic:
    orbit_entries = frontier()
    entry_to_orbit = {
        entry: orbit_index
        for orbit_index, entries in enumerate(orbit_entries)
        for entry in entries
    }
    beta, law = occupation_law()
    weights = dict(law)
    states = tuple(weights)
    state_index = {state: index for index, state in enumerate(states)}
    edges = relevant_edges(orbit_entries, states)
    selected, bases = selected_bases(orbit_entries, states, edges)
    translated = translated_bases(bases, states)
    dimension = len(states) * ORDER
    joint = np.zeros((dimension, dimension))
    by_orbit: list[list[tuple[int, int, float]]] = [[] for _ in orbit_entries]
    activated_entries: set[Entry] = set()
    nonzero_pairs = 0
    matching_pairs = 0

    for left_state, right_state in edges:
        probability_factor = sqrt(weights[left_state] * weights[right_state]) / ORDER
        for left_shift, left in enumerate(translated[left_state]):
            for right_shift, right in enumerate(translated[right_state]):
                supports, entry = difference(left, right)
                orbit_index = entry_to_orbit.get(entry)
                if orbit_index is None:
                    continue
                matching_pairs += 1
                exact = chain_moment(ORDER, supports)
                if not exact:
                    continue
                value = probability_factor * beta ** sum(entry[0]) * float(exact)
                left_index = state_index[left_state] * ORDER + left_shift
                right_index = state_index[right_state] * ORDER + right_shift
                joint[left_index, right_index] += value
                joint[right_index, left_index] += value
                by_orbit[orbit_index].append((left_index, right_index, value))
                activated_entries.add(entry)
                nonzero_pairs += 1

    separate_values = tuple(
        nuclear_from_entries(entries) if entries else 0.0 for entries in by_orbit
    )
    separate = sum(separate_values)
    joint_nuclear = float(np.abs(np.linalg.eigvalsh(joint)).sum())
    current_frontier = sum(float(row["perron_contribution"]) for row in frontier_rows())
    return RowOrbitDiagnostic(
        order=ORDER,
        basis_dimension=dimension,
        occupation_states=len(states),
        translation_orbit_size=ORDER,
        searched_trials=TRIALS,
        selected_trial=selected.trial,
        relevant_occupation_edges=len(edges),
        activated_orbits=sum(bool(entries) for entries in by_orbit),
        activated_entries=len(activated_entries),
        matching_pairs=matching_pairs,
        nonzero_pairs=nonzero_pairs,
        attenuated_separate_nuclear=separate,
        attenuated_joint_nuclear=joint_nuclear,
        attenuated_cancellation_ratio=joint_nuclear / separate if separate else 0.0,
        beta=beta,
        current_frontier_perron_contribution=current_frontier,
        separate_to_current_frontier_ratio=separate / current_frontier,
    )


def artifact_text(result: RowOrbitDiagnostic) -> str:
    payload = {
        "schema": "round4_shared_frontier_row_orbit_v1",
        "result": asdict(result),
        "selection_rule": (
            "maximize exact-nonzero orbits, then exact-nonzero entries, then "
            "absolute exact moment sum, then nonzero pairs"
        ),
        "evidence_label": (
            "exact signed-permutation moments with floating eigendecomposition; "
            "compatible physical lower-witness diagnostic, not theorem bound"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic JSON diagnostic under artifacts/",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "shared_frontier_row_orbit.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "shared frontier q32 row orbit: "
        f"basis={result.basis_dimension},"
        f"selected_trial={result.selected_trial},"
        f"orbits={result.activated_orbits},"
        f"entries={result.activated_entries},"
        f"matching_pairs={result.matching_pairs},"
        f"nonzero_pairs={result.nonzero_pairs},"
        f"separate={result.attenuated_separate_nuclear:.12g},"
        f"joint={result.attenuated_joint_nuclear:.12g},"
        f"ratio={result.attenuated_cancellation_ratio:.12g},"
        f"frontier_ratio={result.separate_to_current_frontier_ratio:.12g},"
        "status=physical_lower_diagnostic_not_theorem"
    )


if __name__ == "__main__":
    main()
