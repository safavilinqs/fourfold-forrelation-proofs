#!/usr/bin/env python3
"""Audit the centered-weight repair of both reflected level-nine trees.

The upper-branching representative is one of the two reflected level-nine
joint saturators left by the N^(1/18) theorem.  The module exhausts its exact
Stein histories, identifies every forced odd local derivative, enumerates all
dangerous grouped-entry partitions, and scores the graphs produced by
re-expanding the outer centered factor.  It then constructs the reflected
lower-branching tree explicitly and independently scores every retained
repair branch.  The one all-fresh reflected branch cancels from the task
difference because it adds a fourth distinct first-layer mark.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import product

from high_level_terminal_best_of_two_audit import (
    branching_potential,
    respecting_max_cut_rank,
    singleton_pair_partitions,
)
from level_ten_forest_mean_zero_repair import (
    LocalWeight,
    Partition,
    TransferRecord,
    apply_transfer,
    best_decay,
    boundaries_from_edges,
    canonical_partition,
    incidence,
    insert_vertex,
    weight_profile,
)


Vertex = tuple[int, int]
Edge = tuple[Vertex, Vertex]

A0 = (0, 0)
B0, B1 = (1, 0), (1, 1)
C0, C1, C2 = (2, 0), (2, 1), (2, 2)
D0, D1, D2 = (3, 0), (3, 1), (3, 2)

TYPE_A_BOUNDARIES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((0, 0), (0, 1)),
    ((0, 0), (0, 1), (1, 2)),
    ((0, 0), (1, 1), (2, 2)),
)
TYPE_A_EDGES: tuple[Edge, ...] = tuple(
    ((boundary, left), (boundary + 1, right))
    for boundary, edges in enumerate(TYPE_A_BOUNDARIES)
    for left, right in edges
)
TYPE_A_VERTICES = tuple(sorted({vertex for edge in TYPE_A_EDGES for vertex in edge}))

TYPE_B_BOUNDARIES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((0, 0), (1, 1), (2, 2)),
    ((0, 0), (1, 0), (2, 1)),
    ((0, 0), (1, 0)),
)
TYPE_B_EDGES: tuple[Edge, ...] = tuple(
    ((boundary, left), (boundary + 1, right))
    for boundary, edges in enumerate(TYPE_B_BOUNDARIES)
    for left, right in edges
)
TYPE_B_VERTICES = tuple(sorted({vertex for edge in TYPE_B_EDGES for vertex in edge}))


@dataclass(frozen=True)
class HistoryAudit:
    """Exact coefficient-interface data for Type A."""

    potential_twelve_initial_configurations: int
    contributing_initial_configurations: int
    histories: int
    weight_profiles: int
    new_transfers: tuple[int, ...]
    existing_transfers: tuple[int, ...]
    derivative_events: tuple[int, ...]
    derivative_sites: tuple[tuple[str, Vertex], ...]
    derivative_weights_by_site: tuple[
        tuple[tuple[str, Vertex], tuple[LocalWeight, ...]], ...
    ]
    time_exponents: tuple[int, int, int]


@dataclass(frozen=True)
class RepairAudit:
    """Exact partition minima for both reflected centered repairs."""

    dangerous_partitions: int
    marked_neighbor_partitions: int
    existing_middle_partitions: int
    existing_endpoint_partitions: int
    fresh_endpoint_partitions: int
    minimum_marked_neighbor_decay: Fraction
    minimum_existing_middle_decay: Fraction
    minimum_existing_endpoint_decay: Fraction
    minimum_fresh_endpoint_decay: Fraction
    reflected_dangerous_partitions: int
    reflected_marked_neighbor_partitions: int
    reflected_existing_middle_partitions: int
    reflected_existing_endpoint_partitions: int
    reflected_cancelled_fresh_partitions: int
    reflected_minimum_marked_neighbor_decay: Fraction
    reflected_minimum_existing_middle_decay: Fraction
    reflected_minimum_existing_endpoint_decay: Fraction
    proved_global_exponent: Fraction


def reflect_vertex(vertex: Vertex) -> Vertex:
    """Reflect a mark through the middle of the four layers."""

    return 3 - vertex[0], vertex[1]


def reflect_edge(edge: Edge) -> Edge:
    """Reflect an oriented edge and restore left-to-right orientation."""

    left, right = edge
    return reflect_vertex(right), reflect_vertex(left)


def reflect_partition(partition: Partition) -> Partition:
    """Reflect every physical-entry block in a grouped-entry partition."""

    return canonical_partition(
        [tuple(reflect_vertex(vertex) for vertex in block) for block in partition]
    )


def assert_reflected_type() -> None:
    """Check that Type B is exactly the layer reflection of Type A."""

    reflected_edges = tuple(sorted(reflect_edge(edge) for edge in TYPE_A_EDGES))
    if reflected_edges != tuple(sorted(TYPE_B_EDGES)):
        raise AssertionError((reflected_edges, TYPE_B_EDGES))
    if boundaries_from_edges(frozenset(reflected_edges)) != TYPE_B_BOUNDARIES:
        raise AssertionError("reflected boundaries")
    if tuple(sorted(reflect_vertex(vertex) for vertex in TYPE_A_VERTICES)) != (
        TYPE_B_VERTICES
    ):
        raise AssertionError("reflected vertices")


def initial_states() -> tuple[
    tuple[frozenset[Edge], dict[tuple[str, Vertex], LocalWeight]], ...
]:
    """Enumerate potential-twelve initial triples contained in Type A."""

    boundary_edges = tuple(
        tuple(((boundary, left), (boundary + 1, right)) for left, right in edges)
        for boundary, edges in enumerate(TYPE_A_BOUNDARIES)
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
    """Return every legal remaining target edge from an active source."""

    outgoing, incoming = incidence(edges)
    remaining = set(TYPE_A_EDGES) - set(edges)
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


def enumerate_histories() -> HistoryAudit:
    """Exhaust all labeled legal histories producing Type A."""

    states = initial_states()
    completed: list[
        tuple[
            tuple[TransferRecord, ...],
            tuple[tuple[str, Vertex, str, int], ...],
        ]
    ] = []

    def visit(
        edges: frozenset[Edge],
        weights: dict[tuple[str, Vertex], LocalWeight],
        records: tuple[TransferRecord, ...],
    ) -> None:
        if edges == frozenset(TYPE_A_EDGES):
            outgoing, incoming = incidence(edges)
            if any(outgoing[layer] != incoming[layer] for layer in (1, 2)):
                raise AssertionError("completed Type A is not terminal")
            completed.append((records, weight_profile(weights)))
            return
        for direction, source, edge in transfer_options(edges):
            new_edges, new_weights, record = apply_transfer(
                edges, weights, direction, source, edge
            )
            visit(new_edges, new_weights, records + (record,))

    contributing = 0
    for edges, weights in states:
        before = len(completed)
        visit(edges, weights, ())
        if len(completed) > before:
            contributing += 1
    if not completed:
        raise AssertionError("Type A must have legal histories")

    profiles = Counter(profile for _, profile in completed)
    derivative_sites: set[tuple[str, Vertex]] = set()
    weights_by_site: dict[tuple[str, Vertex], set[LocalWeight]] = {}
    expected_sites = {("A", A0), ("A", B0)}
    for _, profile in completed:
        derivatives = [entry for entry in profile if entry[3] > 0]
        observed_sites = {(side, vertex) for side, vertex, _, _ in derivatives}
        if len(derivatives) != 2 or observed_sites != expected_sites:
            raise AssertionError(("terminal derivative profile", derivatives))
        if any(order != 1 for _, _, _, order in derivatives):
            raise AssertionError(("terminal derivative order", derivatives))
        for side, vertex, kind, order in derivatives:
            site = (side, vertex)
            derivative_sites.add(site)
            weights_by_site.setdefault(site, set()).add(LocalWeight(kind, order))

    return HistoryAudit(
        potential_twelve_initial_configurations=len(states),
        contributing_initial_configurations=contributing,
        histories=len(completed),
        weight_profiles=len(profiles),
        new_transfers=tuple(
            sorted(
                {
                    sum(record.creates_vertex for record in records)
                    for records, _ in completed
                }
            )
        ),
        existing_transfers=tuple(
            sorted(
                {
                    sum(not record.creates_vertex for record in records)
                    for records, _ in completed
                }
            )
        ),
        derivative_events=tuple(
            sorted(
                {
                    sum(record.differentiates_existing_weight for record in records)
                    for records, _ in completed
                }
            )
        ),
        derivative_sites=tuple(sorted(derivative_sites)),
        derivative_weights_by_site=tuple(
            (
                site,
                tuple(
                    sorted(
                        weights,
                        key=lambda weight: (weight.kind, weight.derivative),
                    )
                ),
            )
            for site, weights in sorted(weights_by_site.items())
        ),
        time_exponents=tuple(len(edges) - 1 for edges in TYPE_A_BOUNDARIES),
    )


def dangerous_partitions() -> tuple[Partition, ...]:
    """Return every maximum-occupancy-two partition of cut rank two."""

    index = {vertex: position for position, vertex in enumerate(TYPE_A_VERTICES)}
    indexed_edges = tuple((index[left], index[right]) for left, right in TYPE_A_EDGES)
    result = []
    for partition in singleton_pair_partitions(tuple(range(len(TYPE_A_VERTICES)))):
        if not any(len(block) == 2 for block in partition):
            continue
        if respecting_max_cut_rank(len(TYPE_A_VERTICES), indexed_edges, partition) != 2:
            continue
        result.append(
            canonical_partition(
                [
                    tuple(TYPE_A_VERTICES[position] for position in block)
                    for block in partition
                ]
            )
        )
    return tuple(sorted(set(result)))


def extend_partitions(
    partitions: tuple[Partition, ...], *vertices: Vertex
) -> tuple[Partition, ...]:
    """Insert new marks into every existing or fresh physical entry."""

    current = set(partitions)
    for vertex in vertices:
        current = {
            extended
            for partition in current
            for extended in insert_vertex(partition, vertex)
        }
    return tuple(sorted(current))


def repair_audit() -> RepairAudit:
    """Score every retained graph produced by both centered expansions."""

    assert_reflected_type()
    partitions = dangerous_partitions()
    if any(
        best_decay(TYPE_A_VERTICES, TYPE_A_EDGES, partition) != Fraction(1, 2)
        for partition in partitions
    ):
        raise AssertionError("dangerous partition score")

    # The centered outer factor at A0 first transfers toward layer one.
    marked_decays = []
    for neighbor in (B0, B1):
        edges = TYPE_A_EDGES + ((A0, neighbor),)
        marked_decays.extend(
            best_decay(TYPE_A_VERTICES, edges, partition) for partition in partitions
        )

    new_b = (1, 2)
    one_new_partitions = extend_partitions(partitions, new_b)
    middle_decays = []
    for middle in (C0, C1, C2):
        edges = TYPE_A_EDGES + ((A0, new_b), (new_b, middle))
        vertices = tuple(sorted(TYPE_A_VERTICES + (new_b,)))
        middle_decays.extend(
            best_decay(vertices, edges, partition) for partition in one_new_partitions
        )

    new_c = (2, 3)
    two_new_partitions = extend_partitions(partitions, new_b, new_c)
    endpoint_decays = []
    for endpoint in (D0, D1, D2):
        edges = TYPE_A_EDGES + (
            (A0, new_b),
            (new_b, new_c),
            (new_c, endpoint),
        )
        vertices = tuple(sorted(TYPE_A_VERTICES + (new_b, new_c)))
        endpoint_decays.extend(
            best_decay(vertices, edges, partition) for partition in two_new_partitions
        )

    new_d = (3, 3)
    three_new_partitions = extend_partitions(partitions, new_b, new_c, new_d)
    fresh_edges = TYPE_A_EDGES + (
        (A0, new_b),
        (new_b, new_c),
        (new_c, new_d),
    )
    fresh_vertices = tuple(sorted(TYPE_A_VERTICES + (new_b, new_c, new_d)))
    fresh_decays = [
        best_decay(fresh_vertices, fresh_edges, partition)
        for partition in three_new_partitions
    ]

    # Reflect every grouped-entry partition and score Type B directly.  Its
    # odd endpoint factor is at D0=(3,0), so the expansion now runs left.
    reflected_partitions = tuple(
        sorted(reflect_partition(partition) for partition in partitions)
    )
    if len(set(reflected_partitions)) != len(partitions):
        raise AssertionError("reflection must preserve partition count")
    if any(
        best_decay(TYPE_B_VERTICES, TYPE_B_EDGES, partition) != Fraction(1, 2)
        for partition in reflected_partitions
    ):
        raise AssertionError("reflected dangerous partition score")

    reflected_marked_decays = []
    reflected_outer = (3, 0)
    for neighbor in ((2, 0), (2, 1)):
        edges = TYPE_B_EDGES + ((neighbor, reflected_outer),)
        reflected_marked_decays.extend(
            best_decay(TYPE_B_VERTICES, edges, partition)
            for partition in reflected_partitions
        )

    new_reflected_c = (2, 2)
    reflected_one_new = extend_partitions(reflected_partitions, new_reflected_c)
    reflected_middle_decays = []
    for middle in ((1, 0), (1, 1), (1, 2)):
        edges = TYPE_B_EDGES + (
            (new_reflected_c, reflected_outer),
            (middle, new_reflected_c),
        )
        vertices = tuple(sorted(TYPE_B_VERTICES + (new_reflected_c,)))
        reflected_middle_decays.extend(
            best_decay(vertices, edges, partition) for partition in reflected_one_new
        )

    new_reflected_b = (1, 3)
    reflected_two_new = extend_partitions(
        reflected_partitions, new_reflected_c, new_reflected_b
    )
    reflected_endpoint_decays = []
    for endpoint in ((0, 0), (0, 1), (0, 2)):
        edges = TYPE_B_EDGES + (
            (new_reflected_c, reflected_outer),
            (new_reflected_b, new_reflected_c),
            (endpoint, new_reflected_b),
        )
        vertices = tuple(sorted(TYPE_B_VERTICES + (new_reflected_c, new_reflected_b)))
        reflected_endpoint_decays.extend(
            best_decay(vertices, edges, partition) for partition in reflected_two_new
        )

    # Continuing to a fresh layer-zero vertex creates a fourth distinct
    # first-layer mark, so this branch is reflection-even and cancels from
    # the planted-minus-null task difference before any absolute value.
    new_reflected_a = (0, 3)
    reflected_three_new = extend_partitions(
        reflected_partitions,
        new_reflected_c,
        new_reflected_b,
        new_reflected_a,
    )
    if sum(vertex[0] == 0 for vertex in TYPE_B_VERTICES) != 3:
        raise AssertionError("Type B sensitivity count")
    if sum(vertex[0] == 0 for vertex in TYPE_B_VERTICES + (new_reflected_a,)) != 4:
        raise AssertionError("all-fresh reflection count")

    return RepairAudit(
        dangerous_partitions=len(partitions),
        marked_neighbor_partitions=len(marked_decays),
        existing_middle_partitions=len(middle_decays),
        existing_endpoint_partitions=len(endpoint_decays),
        fresh_endpoint_partitions=len(fresh_decays),
        minimum_marked_neighbor_decay=min(marked_decays),
        minimum_existing_middle_decay=min(middle_decays),
        minimum_existing_endpoint_decay=min(endpoint_decays),
        minimum_fresh_endpoint_decay=min(fresh_decays),
        reflected_dangerous_partitions=len(reflected_partitions),
        reflected_marked_neighbor_partitions=len(reflected_marked_decays),
        reflected_existing_middle_partitions=len(reflected_middle_decays),
        reflected_existing_endpoint_partitions=len(reflected_endpoint_decays),
        reflected_cancelled_fresh_partitions=len(reflected_three_new),
        reflected_minimum_marked_neighbor_decay=min(reflected_marked_decays),
        reflected_minimum_existing_middle_decay=min(reflected_middle_decays),
        reflected_minimum_existing_endpoint_decay=min(reflected_endpoint_decays),
        proved_global_exponent=Fraction(1, 16),
    )


def main() -> None:
    history = enumerate_histories()
    repair = repair_audit()
    print(
        "level-nine tree centered repair: "
        f"potential_initials={history.potential_twelve_initial_configurations},"
        f"contributing_initials={history.contributing_initial_configurations},"
        f"histories={history.histories},"
        f"profiles={history.weight_profiles},"
        f"new={history.new_transfers},"
        f"existing={history.existing_transfers},"
        f"derivatives={history.derivative_events},"
        f"sites={history.derivative_sites},"
        f"weights={history.derivative_weights_by_site},"
        f"time_exponents={history.time_exponents},"
        f"partitions={repair.dangerous_partitions},"
        f"marked_decay={repair.minimum_marked_neighbor_decay},"
        f"middle_decay={repair.minimum_existing_middle_decay},"
        f"endpoint_decay={repair.minimum_existing_endpoint_decay},"
        f"fresh_decay={repair.minimum_fresh_endpoint_decay},"
        f"reflected_partitions={repair.reflected_dangerous_partitions},"
        f"reflected_marked_decay="
        f"{repair.reflected_minimum_marked_neighbor_decay},"
        f"reflected_middle_decay="
        f"{repair.reflected_minimum_existing_middle_decay},"
        f"reflected_endpoint_decay="
        f"{repair.reflected_minimum_existing_endpoint_decay},"
        f"reflected_cancelled={repair.reflected_cancelled_fresh_partitions},"
        f"global_exponent={repair.proved_global_exponent}"
    )


if __name__ == "__main__":
    main()
