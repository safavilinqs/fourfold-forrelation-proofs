#!/usr/bin/env python3
"""Audit the centered-weight repair of the unique level-ten saturator.

The dangerous terminal graph has a six-vertex tree and a four-layer path.
Every legal Stein history producing this graph contains exactly one
same-incidence transfer into the degree-three middle vertex.  Consequently
its local analytic weight contains one odd centered derivative, either
``gamma' = psi''`` or the derivative of the even Stein kernel.

Because ``psi''`` is odd and Gaussian-centered, one more Stein identity
adds an edge from that middle vertex.  If the neighbor is already marked,
the resulting graph is checked directly below.  If the neighbor is new, a
second transfer either creates a new first-layer mark (reflection-even and
therefore absent from the task difference) or hits an existing first-layer
mark.  The latter graphs are also checked below for every dangerous
physical-entry partition.

This file certifies only finite combinatorics and exact GF(2) cut ranks.  The
companion note states the analytic centered-Stein identity and combines it
with the existing grouped-entry projective contraction.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import product

from high_level_terminal_best_of_two_audit import (
    branching_potential,
    cut_rank,
)


Vertex = tuple[int, int]
Edge = tuple[Vertex, Vertex]
Block = tuple[Vertex, ...]
Partition = tuple[Block, ...]

A0, A1, A2 = (0, 0), (0, 1), (0, 2)
B0, B1, B2 = (1, 0), (1, 1), (1, 2)
C0, C1 = (2, 0), (2, 1)
D0, D1 = (3, 0), (3, 1)

TYPE_C_BOUNDARIES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((0, 0), (1, 1), (2, 2)),
    ((0, 0), (1, 0), (2, 1)),
    ((0, 0), (1, 1)),
)

TYPE_C_EDGES: tuple[Edge, ...] = tuple(
    ((boundary, left), (boundary + 1, right))
    for boundary, edges in enumerate(TYPE_C_BOUNDARIES)
    for left, right in edges
)
TYPE_C_VERTICES = tuple(sorted({vertex for edge in TYPE_C_EDGES for vertex in edge}))
STRONG_COMPONENT = (A0, A1, B0, B1, C0, D0)
PATH_COMPONENT = (A2, B2, C1, D1)


@dataclass(frozen=True)
class LocalWeight:
    """One local gamma or Stein-kernel factor and its derivative order."""

    kind: str
    derivative: int = 0


@dataclass(frozen=True)
class TransferRecord:
    """One legal edge-adding transfer along a fixed labeled history."""

    edge: Edge
    direction: str
    source: Vertex
    neighbor: Vertex
    creates_vertex: bool
    differentiates_existing_weight: bool


@dataclass(frozen=True)
class HistoryAudit:
    """Exact coefficient-interface summary for all Type-C histories."""

    potential_twelve_initial_configurations: int
    contributing_initial_configurations: int
    histories: int
    weight_profiles: int
    new_transfers: tuple[int, ...]
    existing_transfers: tuple[int, ...]
    derivative_events: tuple[int, ...]
    common_derivative_site: tuple[str, Vertex]
    common_derivative_weights: tuple[LocalWeight, ...]
    time_exponents: tuple[int, int, int]


@dataclass(frozen=True)
class RepairAudit:
    """Physical partition and graph-decay verdict for the repair branches."""

    dangerous_partitions: int
    duplicate_branches: int
    bridge_partitions: int
    new_neighbor_partitions: int
    minimum_duplicate_decay: Fraction
    minimum_bridge_decay: Fraction
    minimum_existing_outer_decay: Fraction
    repaired_level_ten_decay: Fraction
    repaired_level_eleven_decay: Fraction
    proved_global_exponent: Fraction


def boundaries_from_edges(
    edges: frozenset[Edge],
) -> tuple[tuple[tuple[int, int], ...], ...]:
    """Return the set-valued boundary representation used by the potential."""

    result = []
    for boundary in range(3):
        result.append(
            tuple(
                sorted(
                    (left[1], right[1])
                    for left, right in edges
                    if left[0] == boundary and right[0] == boundary + 1
                )
            )
        )
    return tuple(result)


def incidence(edges: frozenset[Edge]) -> tuple[list[set[int]], list[set[int]]]:
    """Return outgoing and incoming marked-coordinate sets."""

    outgoing = [set() for _ in range(4)]
    incoming = [set() for _ in range(4)]
    for left, right in edges:
        outgoing[left[0]].add(left[1])
        incoming[right[0]].add(right[1])
    return outgoing, incoming


def initial_states() -> tuple[
    tuple[frozenset[Edge], dict[tuple[str, Vertex], LocalWeight]], ...
]:
    """Enumerate the high-level initial edge choices contained in Type C."""

    boundary_edges = tuple(
        tuple(((boundary, left), (boundary + 1, right)) for left, right in edges)
        for boundary, edges in enumerate(TYPE_C_BOUNDARIES)
    )
    states = []
    for chosen in product(*boundary_edges):
        edge_set = frozenset(chosen)
        if branching_potential(boundaries_from_edges(edge_set)) != 12:
            continue
        weights: dict[tuple[str, Vertex], LocalWeight] = {}
        for left, right in chosen:
            weights[("A", left)] = LocalWeight("gamma")
            weights[("B", right)] = LocalWeight("gamma")
        states.append((edge_set, weights))
    return tuple(states)


def transfer_options(edges: frozenset[Edge]) -> tuple[tuple[str, Vertex, Edge], ...]:
    """Return legal remaining Type-C edges from every active source."""

    outgoing, incoming = incidence(edges)
    remaining = set(TYPE_C_EDGES) - set(edges)
    result = []
    for layer in (1, 2):
        for coordinate in sorted(incoming[layer] - outgoing[layer]):
            source = (layer, coordinate)
            result.extend(
                ("right", source, edge) for edge in remaining if edge[0] == source
            )
        for coordinate in sorted(outgoing[layer] - incoming[layer]):
            source = (layer, coordinate)
            result.extend(
                ("left", source, edge) for edge in remaining if edge[1] == source
            )
    return tuple(sorted(result))


def apply_transfer(
    edges: frozenset[Edge],
    weights: dict[tuple[str, Vertex], LocalWeight],
    direction: str,
    source: Vertex,
    edge: Edge,
) -> tuple[frozenset[Edge], dict[tuple[str, Vertex], LocalWeight], TransferRecord]:
    """Apply one transfer and update the exact local function labels."""

    outgoing, incoming = incidence(edges)
    occupied = {vertex for old_edge in edges for vertex in old_edge}
    updated = dict(weights)
    if direction == "right":
        if source != edge[0]:
            raise AssertionError((direction, source, edge))
        neighbor = edge[1]
        updated[("A", source)] = LocalWeight("stein")
        key = ("B", neighbor)
        differentiates = neighbor[1] in incoming[neighbor[0]]
    else:
        if source != edge[1]:
            raise AssertionError((direction, source, edge))
        neighbor = edge[0]
        updated[("B", source)] = LocalWeight("stein")
        key = ("A", neighbor)
        differentiates = neighbor[1] in outgoing[neighbor[0]]

    if differentiates:
        old = updated[key]
        updated[key] = LocalWeight(old.kind, old.derivative + 1)
    else:
        if key in updated:
            raise AssertionError(("unexpected existing function", key))
        updated[key] = LocalWeight("gamma")

    new_edges = edges | {edge}
    before = branching_potential(boundaries_from_edges(edges))
    after = branching_potential(boundaries_from_edges(new_edges))
    creates = neighbor not in occupied
    if creates and before != after:
        raise AssertionError(("fresh potential", before, after, edge))
    if not creates and not after < before:
        raise AssertionError(("existing potential", before, after, edge))
    return (
        new_edges,
        updated,
        TransferRecord(edge, direction, source, neighbor, creates, differentiates),
    )


def weight_profile(
    weights: dict[tuple[str, Vertex], LocalWeight],
) -> tuple[tuple[str, Vertex, str, int], ...]:
    """Canonicalize a terminal local-weight profile."""

    return tuple(
        sorted(
            (side, vertex, weight.kind, weight.derivative)
            for (side, vertex), weight in weights.items()
        )
    )


def enumerate_histories() -> HistoryAudit:
    """Exhaust every labeled initial choice and legal transfer order."""

    states = initial_states()
    completed: list[
        tuple[tuple[TransferRecord, ...], tuple[tuple[str, Vertex, str, int], ...]]
    ] = []

    def visit(
        edges: frozenset[Edge],
        weights: dict[tuple[str, Vertex], LocalWeight],
        records: tuple[TransferRecord, ...],
    ) -> None:
        if edges == frozenset(TYPE_C_EDGES):
            outgoing, incoming = incidence(edges)
            if any(outgoing[layer] != incoming[layer] for layer in (1, 2)):
                raise AssertionError("completed Type C is not terminal")
            completed.append((records, weight_profile(weights)))
            return
        options = transfer_options(edges)
        if not options:
            return
        for direction, source, edge in options:
            new_edges, new_weights, record = apply_transfer(
                edges, weights, direction, source, edge
            )
            visit(new_edges, new_weights, records + (record,))

    contributing_initials = 0
    for edges, weights in states:
        previous_histories = len(completed)
        visit(edges, weights, ())
        if len(completed) > previous_histories:
            contributing_initials += 1

    if not completed:
        raise AssertionError("Type C must have legal histories")
    profiles = Counter(profile for _, profile in completed)
    new_counts = tuple(
        sorted(
            {
                sum(record.creates_vertex for record in records)
                for records, _ in completed
            }
        )
    )
    existing_counts = tuple(
        sorted(
            {
                sum(not record.creates_vertex for record in records)
                for records, _ in completed
            }
        )
    )
    derivative_counts = tuple(
        sorted(
            {
                sum(record.differentiates_existing_weight for record in records)
                for records, _ in completed
            }
        )
    )
    derivative_sites = set()
    derivative_weights = set()
    for _, profile in completed:
        derivatives = [entry for entry in profile if entry[3] > 0]
        if len(derivatives) != 1:
            raise AssertionError(("terminal derivative profile", derivatives))
        side, vertex, kind, order = derivatives[0]
        derivative_sites.add((side, vertex))
        derivative_weights.add(LocalWeight(kind, order))
    if len(derivative_sites) != 1:
        raise AssertionError((derivative_sites, derivative_weights))

    edge_counts = tuple(len(edges) for edges in TYPE_C_BOUNDARIES)
    time_exponents = tuple(count - 1 for count in edge_counts)
    return HistoryAudit(
        potential_twelve_initial_configurations=len(states),
        contributing_initial_configurations=contributing_initials,
        histories=len(completed),
        weight_profiles=len(profiles),
        new_transfers=new_counts,
        existing_transfers=existing_counts,
        derivative_events=derivative_counts,
        common_derivative_site=next(iter(derivative_sites)),
        common_derivative_weights=tuple(
            sorted(
                derivative_weights, key=lambda weight: (weight.kind, weight.derivative)
            )
        ),
        time_exponents=time_exponents,
    )


def canonical_partition(blocks: list[list[Vertex]] | tuple[Block, ...]) -> Partition:
    """Return a deterministic tuple representation of physical entries."""

    normalized = [tuple(sorted(block)) for block in blocks if block]
    return tuple(sorted(normalized, key=lambda block: (block[0], len(block), block)))


def dangerous_partitions() -> tuple[Partition, ...]:
    """Enumerate every global partition inducing the Type-C worst score.

    The strong component has exactly two rank-one induced partitions.  Each
    path vertex is singleton within its own component, but may share an
    entry with one strong-component block.
    """

    strong_partitions = (
        ((A0, B0), (A1, B1), (C0,), (D0,)),
        ((A0, B0), (A1, B1), (C0, D0)),
    )
    results: set[Partition] = set()
    path_vertices = PATH_COMPONENT
    for strong in strong_partitions:
        blocks = [list(block) for block in strong]

        def assign(index: int, available: tuple[int, ...]) -> None:
            if index == len(path_vertices):
                results.add(canonical_partition(blocks))
                return
            vertex = path_vertices[index]
            blocks.append([vertex])
            assign(index + 1, available)
            blocks.pop()
            for position, block_index in enumerate(available):
                blocks[block_index].append(vertex)
                assign(index + 1, available[:position] + available[position + 1 :])
                blocks[block_index].pop()

        assign(0, tuple(range(len(strong))))
    return tuple(sorted(results))


def graph_components(
    vertices: tuple[Vertex, ...], edges: tuple[Edge, ...]
) -> tuple[tuple[Vertex, ...], ...]:
    """Return ordinary connected components, retaining parallel edges."""

    adjacency: dict[Vertex, set[Vertex]] = defaultdict(set)
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    remaining = set(vertices)
    components = []
    while remaining:
        start = min(remaining)
        stack = [start]
        component = {start}
        remaining.remove(start)
        while stack:
            vertex = stack.pop()
            for neighbor in adjacency[vertex]:
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    component.add(neighbor)
                    stack.append(neighbor)
        components.append(tuple(sorted(component)))
    return tuple(components)


def induced_partition(component: tuple[Vertex, ...], partition: Partition) -> Partition:
    """Restrict global physical entries to one graph component."""

    component_set = set(component)
    return canonical_partition(
        [[vertex for vertex in block if vertex in component_set] for block in partition]
    )


def component_cut_rank(
    component: tuple[Vertex, ...],
    edges: tuple[Edge, ...],
    partition: Partition,
) -> int:
    """Maximum GF(2) adjacency rank over entry-respecting cuts."""

    index = {vertex: position for position, vertex in enumerate(component)}
    indexed_edges = tuple(
        (index[left], index[right])
        for left, right in edges
        if left in index and right in index
    )
    blocks = tuple(tuple(index[vertex] for vertex in block) for block in partition)
    ranks = []
    for mask in range(1, (1 << len(blocks)) - 1):
        selected = frozenset(
            vertex
            for block_index, block in enumerate(blocks)
            if mask >> block_index & 1
            for vertex in block
        )
        ranks.append(cut_rank(len(component), indexed_edges, selected))
    return max(ranks)


def best_decay(
    vertices: tuple[Vertex, ...],
    edges: tuple[Edge, ...],
    partition: Partition,
) -> Fraction:
    """Return the better complete assigned and all-projective exponents."""

    assigned_sigma = 0
    projective_decay = Fraction()
    for component in graph_components(vertices, edges):
        component_set = set(component)
        edge_count = sum(
            left in component_set and right in component_set for left, right in edges
        )
        surplus = edge_count - len(component) + 1
        induced = induced_partition(component, partition)
        maximum_occupancy = max(map(len, induced))
        assigned_sigma += surplus + maximum_occupancy - 1
        rank = component_cut_rank(component, edges, induced)
        projective_decay += Fraction(surplus + rank - 1, 2)
    return max(Fraction(assigned_sigma, 2), projective_decay)


def insert_vertex(partition: Partition, vertex: Vertex) -> tuple[Partition, ...]:
    """Place a new mark in any existing physical entry or a fresh entry."""

    results = {canonical_partition([*partition, (vertex,)])}
    for index in range(len(partition)):
        blocks = [list(block) for block in partition]
        blocks[index].append(vertex)
        results.add(canonical_partition(blocks))
    return tuple(sorted(results))


def verify_original_partitions(partitions: tuple[Partition, ...]) -> None:
    """Check that the enumerated partitions are exactly joint saturators."""

    for partition in partitions:
        if best_decay(TYPE_C_VERTICES, TYPE_C_EDGES, partition) != Fraction(1, 2):
            raise AssertionError(("not a Type-C saturator", partition))


def repair_audit() -> RepairAudit:
    """Audit every reflection-sensitive branch after centering psi''."""

    partitions = dangerous_partitions()
    verify_original_partitions(partitions)

    duplicate_decays = []
    for neighbor in (B0, B1):
        duplicate_edges = TYPE_C_EDGES + ((neighbor, C0),)
        duplicate_decays.extend(
            best_decay(TYPE_C_VERTICES, duplicate_edges, partition)
            for partition in partitions
        )

    bridge_edges = TYPE_C_EDGES + ((B2, C0),)
    bridge_decays = [
        best_decay(TYPE_C_VERTICES, bridge_edges, partition) for partition in partitions
    ]

    new_neighbor = (1, 3)
    outer_decays = []
    extended_partitions = 0
    for outer in (A0, A1, A2):
        extension_edges = TYPE_C_EDGES + ((new_neighbor, C0), (outer, new_neighbor))
        extension_vertices = tuple(sorted(TYPE_C_VERTICES + (new_neighbor,)))
        for partition in partitions:
            for extended in insert_vertex(partition, new_neighbor):
                extended_partitions += 1
                outer_decays.append(
                    best_decay(extension_vertices, extension_edges, extended)
                )

    minimum_duplicate = min(duplicate_decays)
    minimum_bridge = min(bridge_decays)
    minimum_outer = min(outer_decays)
    repaired_level_ten = min(minimum_duplicate, minimum_bridge)
    repaired_level_eleven = minimum_outer

    # After Type C is repaired, the two level-nine trees with N^{-1/2}
    # are limiting.  All levels <=8 retain N^{-1/2}; all other levels have
    # their previously audited stronger rows.
    exponent_candidates = [Fraction(1, 2 * level) for level in range(4, 10)]
    exponent_candidates.extend(
        (
            repaired_level_ten / 10,
            repaired_level_eleven / 11,
            Fraction(1, 12),
        )
    )
    return RepairAudit(
        dangerous_partitions=len(partitions),
        duplicate_branches=2 * len(partitions),
        bridge_partitions=len(partitions),
        new_neighbor_partitions=extended_partitions,
        minimum_duplicate_decay=minimum_duplicate,
        minimum_bridge_decay=minimum_bridge,
        minimum_existing_outer_decay=minimum_outer,
        repaired_level_ten_decay=repaired_level_ten,
        repaired_level_eleven_decay=repaired_level_eleven,
        proved_global_exponent=min(exponent_candidates),
    )


def main() -> None:
    history = enumerate_histories()
    repair = repair_audit()
    print(
        "level-ten forest mean-zero repair: "
        f"potential_initials={history.potential_twelve_initial_configurations},"
        f"contributing_initials={history.contributing_initial_configurations},"
        f"histories={history.histories},"
        f"profiles={history.weight_profiles},"
        f"new={history.new_transfers},"
        f"existing={history.existing_transfers},"
        f"derivatives={history.derivative_events},"
        f"derivative_site={history.common_derivative_site},"
        f"derivative_weights={history.common_derivative_weights},"
        f"time_exponents={history.time_exponents},"
        f"partitions={repair.dangerous_partitions},"
        f"duplicate_decay={repair.minimum_duplicate_decay},"
        f"bridge_decay={repair.minimum_bridge_decay},"
        f"outer_decay={repair.minimum_existing_outer_decay},"
        f"global_exponent={repair.proved_global_exponent}"
    )


if __name__ == "__main__":
    main()
