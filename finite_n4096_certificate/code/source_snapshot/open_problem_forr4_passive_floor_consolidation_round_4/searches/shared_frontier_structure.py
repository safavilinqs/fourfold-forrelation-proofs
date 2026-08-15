#!/usr/bin/env python3
"""Factor the leading unresolved frontier into occupation-kernel families.

This is a structural audit, not a coefficient theorem.  It starts from the
generated Round 4 unresolved-orbit ledger, takes the first 51 orbits carrying
at least 90 percent of the unresolved Perron contribution, and enumerates the
exact dose-six occupation pairs on which those entries act.

For a balanced profile/split pair ``(a, s)`` with ``|a| = 2r`` and ``|s| = r``,
every compatible pair is

    n = s + k,       m = (a-s) + k,

where ``k`` is the shared unmarked occupation.  At hard dose six, the degree
ten frontier therefore has ``|k| <= 1`` and the degree twelve frontier has
``k = 0``.  This collapses the nominal 198 entries to a small shared graph.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from json import dumps, loads
from pathlib import Path
import sys
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent / "open_problem_forr4_passive_floor_consolidation_round_3" / "searches"
)
sys.path.insert(0, str(ROUND3_SEARCHES))

from occupation_compatible_sector_optimization import (  # noqa: E402
    multiplicity,
    paired_state,
)
from double_endpoint_occupation_optimization import (  # noqa: E402
    occupation_states,
)


FRONTIER_ORBITS = 51
DOSE = 6

Profile = tuple[int, ...]
Split = tuple[int, ...]
State = tuple[int, ...]
Entry = tuple[Profile, Split]
Edge = tuple[State, State]


def _counter_payload(counter: Counter[int]) -> dict[str, int]:
    return {str(key): counter[key] for key in sorted(counter)}


def frontier_rows() -> tuple[dict[str, str], ...]:
    path = ROOT / "artifacts" / "n1024_unresolved_orbits.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = tuple(csv.DictReader(handle))
    if len(rows) < FRONTIER_ORBITS:
        raise AssertionError(("unresolved orbit rows", len(rows)))
    return rows[:FRONTIER_ORBITS]


def parsed_entries(row: dict[str, str]) -> tuple[Entry, ...]:
    return tuple(
        (
            tuple(int(value) for value in profile),
            tuple(int(value) for value in split),
        )
        for profile, split in loads(row["entries"])
    )


def compatible_edges(entry: Entry, states: tuple[State, ...]) -> tuple[Edge, ...]:
    profile, split = entry
    if sum(profile) != 2 * sum(split):
        raise AssertionError(("entry is not balanced", entry))
    complement = tuple(
        degree - selected for degree, selected in zip(profile, split, strict=True)
    )
    state_set = set(states)
    result: list[Edge] = []
    for state in states:
        if any(
            occupation < selected
            for occupation, selected in zip(state, split, strict=True)
        ):
            continue
        partner = paired_state(state, profile, split)
        if partner not in state_set:
            continue
        if not multiplicity(state, split):
            continue
        if not multiplicity(partner, complement):
            continue
        intersection = tuple(
            occupation - selected
            for occupation, selected in zip(state, split, strict=True)
        )
        partner_intersection = tuple(
            occupation - selected
            for occupation, selected in zip(partner, complement, strict=True)
        )
        if intersection != partner_intersection:
            raise AssertionError((entry, state, partner))
        result.append(tuple(sorted((state, partner))))
    return tuple(result)


def connected_components(edges: Iterable[Edge]) -> tuple[tuple[State, ...], ...]:
    adjacency: dict[State, set[State]] = defaultdict(set)
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    seen: set[State] = set()
    result: list[tuple[State, ...]] = []
    for start in sorted(adjacency):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        component: list[State] = []
        while stack:
            state = stack.pop()
            component.append(state)
            for neighbor in adjacency[state]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        result.append(tuple(sorted(component)))
    return tuple(sorted(result, key=lambda item: (sum(item[0]), -len(item), item)))


def audit() -> dict[str, object]:
    rows = frontier_rows()
    states = tuple(occupation_states())
    all_entries: set[Entry] = set()
    all_edges: set[Edge] = set()
    edge_entries: dict[Edge, set[Entry]] = defaultdict(set)
    edge_orbits: dict[Edge, set[int]] = defaultdict(set)
    edge_patterns: dict[Edge, set[tuple[int, ...]]] = defaultdict(set)
    degree_entries: Counter[int] = Counter()
    degree_terms: Counter[tuple[int, int]] = Counter()
    degree_edges: dict[int, set[Edge]] = defaultdict(set)
    entry_edge_counts: Counter[int] = Counter()
    orbit_edge_counts: Counter[int] = Counter()
    pattern_stats: dict[tuple[int, ...], dict[str, object]] = defaultdict(
        lambda: {
            "orbits": 0,
            "entries": 0,
            "perron_contribution": 0.0,
            "edges": set(),
            "states": set(),
        }
    )

    for orbit_index, row in enumerate(rows):
        entries = parsed_entries(row)
        if not entries:
            raise AssertionError(("empty orbit", orbit_index))
        pattern = tuple(sorted(entries[0][0], reverse=True))
        if any(
            tuple(sorted(profile, reverse=True)) != pattern for profile, _ in entries
        ):
            raise AssertionError(("mixed profile pattern", orbit_index))
        orbit_edges: set[Edge] = set()
        stats = pattern_stats[pattern]
        stats["orbits"] = int(stats["orbits"]) + 1
        stats["entries"] = int(stats["entries"]) + len(entries)
        stats["perron_contribution"] = float(stats["perron_contribution"]) + float(
            row["perron_contribution"]
        )

        for entry in entries:
            profile, split = entry
            degree = sum(profile)
            all_entries.add(entry)
            degree_entries[degree] += 1
            edges = compatible_edges(entry, states)
            entry_edge_counts[len(edges)] += 1
            for edge in edges:
                left, _ = edge
                intersection_size = sum(left) - sum(split)
                degree_terms[(degree, intersection_size)] += 1
                all_edges.add(edge)
                orbit_edges.add(edge)
                degree_edges[degree].add(edge)
                edge_entries[edge].add(entry)
                edge_orbits[edge].add(orbit_index)
                edge_patterns[edge].add(pattern)
                cast_edges = stats["edges"]
                cast_states = stats["states"]
                if not isinstance(cast_edges, set) or not isinstance(cast_states, set):
                    raise AssertionError("internal pattern accumulator type")
                cast_edges.add(edge)
                cast_states.update(edge)
        orbit_edge_counts[len(orbit_edges)] += 1

    summary = loads(
        (ROOT / "artifacts" / "finite_size_ledger_summary.json").read_text(
            encoding="utf-8"
        )
    )
    total_unresolved = float(
        summary["impact_frontier"]["total_unresolved_perron_contribution"]
    )
    frontier_contribution = sum(float(row["perron_contribution"]) for row in rows)
    components = connected_components(all_edges)
    patterns = []
    for pattern, stats in sorted(
        pattern_stats.items(),
        key=lambda item: (-float(item[1]["perron_contribution"]), item[0]),
    ):
        edges = stats["edges"]
        pattern_states = stats["states"]
        if not isinstance(edges, set) or not isinstance(pattern_states, set):
            raise AssertionError("internal pattern payload type")
        contribution = float(stats["perron_contribution"])
        patterns.append(
            {
                "sorted_profile": list(pattern),
                "degree": sum(pattern),
                "orbits": int(stats["orbits"]),
                "entries": int(stats["entries"]),
                "occupation_edges": len(edges),
                "occupation_states": len(pattern_states),
                "perron_contribution": contribution,
                "fraction_of_frontier_contribution": contribution
                / frontier_contribution,
            }
        )

    degree_payload = []
    for degree in sorted(degree_entries):
        intersections = {
            intersection: degree_terms[(degree, intersection)]
            for candidate_degree, intersection in degree_terms
            if candidate_degree == degree
        }
        degree_payload.append(
            {
                "degree": degree,
                "balanced_split_degree": degree // 2,
                "entries": degree_entries[degree],
                "compatible_terms": sum(intersections.values()),
                "shared_intersection_term_counts": {
                    str(key): intersections[key] for key in sorted(intersections)
                },
                "occupation_layers": sorted(
                    {degree // 2 + key for key in intersections}
                ),
                "unique_occupation_edges": len(degree_edges[degree]),
            }
        )

    component_payload = [
        {
            "occupation_layer": sum(component[0]),
            "states": len(component),
        }
        for component in components
    ]
    return {
        "schema": "round4_shared_frontier_structure_v1",
        "evidence_label": (
            "exact combinatorial occupation-kernel audit; routing structure, "
            "not an arbitrary-law coefficient bound"
        ),
        "dose": DOSE,
        "frontier": {
            "orbits": len(rows),
            "entries": len(all_entries),
            "perron_contribution": frontier_contribution,
            "fraction_of_unresolved_perron_contribution": frontier_contribution
            / total_unresolved,
            "unordered_profile_patterns": len(patterns),
        },
        "degree_layers": degree_payload,
        "kernel": {
            "occupation_states": len({state for edge in all_edges for state in edge}),
            "occupation_edges": len(all_edges),
            "connected_components": len(components),
            "components": component_payload,
            "cross_degree_shared_edges": len(degree_edges[10] & degree_edges[12]),
            "entries_per_edge": _counter_payload(
                Counter(len(entries) for entries in edge_entries.values())
            ),
            "orbits_per_edge": _counter_payload(
                Counter(len(orbits) for orbits in edge_orbits.values())
            ),
            "patterns_per_edge": _counter_payload(
                Counter(len(patterns_) for patterns_ in edge_patterns.values())
            ),
            "edges_per_entry": _counter_payload(entry_edge_counts),
            "edges_per_orbit": _counter_payload(orbit_edge_counts),
        },
        "profile_patterns": patterns,
    }


def artifact_text(result: dict[str, object]) -> str:
    return dumps(result, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic JSON audit under artifacts/",
    )
    arguments = parser.parse_args()
    result = audit()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "shared_frontier_structure.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    frontier = result["frontier"]
    kernel = result["kernel"]
    if not isinstance(frontier, dict) or not isinstance(kernel, dict):
        raise AssertionError("malformed audit payload")
    print(
        "shared frontier structure: "
        f"orbits={frontier['orbits']},"
        f"entries={frontier['entries']},"
        f"fraction={frontier['fraction_of_unresolved_perron_contribution']:.12g},"
        f"patterns={frontier['unordered_profile_patterns']},"
        f"states={kernel['occupation_states']},"
        f"edges={kernel['occupation_edges']},"
        f"components={kernel['connected_components']},"
        "status=exact_structure_not_coefficient_theorem"
    )


if __name__ == "__main__":
    main()
