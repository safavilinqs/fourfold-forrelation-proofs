#!/usr/bin/env python3
"""Audit the full terminal image and repair its low-level saturators.

The high-level enumeration starts from the potential-twelve initial state.
For levels at most eight one must also include the two potential-eight
initial collision patterns and the potential-four path.  This module
enumerates the union of all four initial states, applies the same safe
assigned/projective audit, and then resolves the unique limiting types at
levels seven and six by their forced odd centered local derivatives.

Together with the accepted level-nine and level-ten centered repairs, the
result improves the global passive exponent from 1/16 to 1/12.  The final
limiting row is the level-twelve N^-1 contraction.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import product

from high_level_terminal_best_of_two_audit import (
    Boundaries,
    TerminalState,
    branching_potential,
    canonicalize,
    respecting_max_cut_rank,
    sigma_one_audit,
    singleton_pair_partitions,
    terminal_state,
    transfer_children,
)
from level_nine_tree_centered_repair import (
    extend_partitions,
    repair_audit as level_nine_repair_audit,
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
    repair_audit as level_ten_repair_audit,
    weight_profile,
)


Vertex = tuple[int, int]
Edge = tuple[Vertex, Vertex]

# The four canonical initial covariance patterns: no internal collision,
# a collision in either internal layer, and collisions in both layers.
INITIAL_BOUNDARIES: tuple[Boundaries, ...] = (
    (((0, 0),), ((1, 0),), ((1, 0),)),
    (((0, 0),), ((0, 0),), ((1, 0),)),
    (((0, 0),), ((1, 0),), ((0, 0),)),
    (((0, 0),), ((0, 0),), ((0, 0),)),
)

# The only safe-norm saturator at level seven after all initial collision
# patterns are included: two four-layer paths sharing their first vertex.
LEVEL_SEVEN_BOUNDARIES: Boundaries = (
    ((0, 0), (0, 1)),
    ((0, 0), (1, 1)),
    ((0, 0), (1, 1)),
)

# The only safe-norm saturator at level six.  Its centered factor is internal
# and its fresh outer branch cancels by first-layer reflection.
LEVEL_SIX_BOUNDARIES: Boundaries = (
    ((0, 0),),
    ((0, 0), (0, 1)),
    ((0, 0), (1, 1)),
)


@dataclass(frozen=True)
class FullImageAudit:
    """Exact terminal-image and pre-repair norm summary."""

    initial_potentials: tuple[int, ...]
    reachable_states: int
    terminal_types: int
    sensitive_terminal_types: int
    terminal_types_by_level: tuple[tuple[int, int], ...]
    sensitive_types_by_level: tuple[tuple[int, int], ...]
    minimum_safe_decay_by_level: tuple[tuple[int, Fraction], ...]
    sensitive_level_eight_types: int
    level_eight_minimum_decay: Fraction
    level_seven_saturators: int
    level_six_saturators: int


@dataclass(frozen=True)
class CoefficientAudit:
    """All exact coefficient histories for one low-level target."""

    level: int
    potential_initial_configurations: tuple[tuple[int, int], ...]
    contributing_initial_configurations: tuple[tuple[int, int], ...]
    histories_by_initial_potential: tuple[tuple[int, int], ...]
    histories: int
    weight_profiles: int
    fresh_transfers_by_initial_potential: tuple[tuple[int, tuple[int, ...]], ...]
    existing_transfers_by_initial_potential: tuple[tuple[int, tuple[int, ...]], ...]
    derivative_events_by_initial_potential: tuple[tuple[int, tuple[int, ...]], ...]
    common_derivative_sites: tuple[tuple[str, Vertex], ...]
    common_derivative_weights: tuple[
        tuple[tuple[str, Vertex], tuple[LocalWeight, ...]], ...
    ]
    time_exponents: tuple[int, int, int]


@dataclass(frozen=True)
class LowLevelRepairAudit:
    """Exact physical-partition minima for the level-seven/six repairs."""

    level_seven_dangerous_partitions: int
    level_seven_marked_neighbor_cases: int
    level_seven_existing_middle_cases: int
    level_seven_existing_endpoint_cases: int
    level_seven_fresh_endpoint_cases: int
    level_seven_minimum_marked_neighbor_decay: Fraction
    level_seven_minimum_existing_middle_decay: Fraction
    level_seven_minimum_existing_endpoint_decay: Fraction
    level_seven_minimum_fresh_endpoint_decay: Fraction
    level_six_dangerous_partitions: int
    level_six_marked_outer_cases: int
    level_six_cancelled_fresh_outer_cases: int
    level_six_minimum_marked_outer_decay: Fraction
    proved_global_exponent: Fraction


@dataclass(frozen=True)
class CompleteAudit:
    """Combined full-image, coefficient, repair, and exponent certificate."""

    image: FullImageAudit
    level_seven_coefficients: CoefficientAudit
    level_six_coefficients: CoefficientAudit
    repair: LowLevelRepairAudit


def edges_from_boundaries(boundaries: Boundaries) -> tuple[Edge, ...]:
    """Return ordinary oriented graph edges."""

    return tuple(
        ((boundary, left), (boundary + 1, right))
        for boundary, boundary_edges in enumerate(boundaries)
        for left, right in boundary_edges
    )


def vertices_from_edges(edges: tuple[Edge, ...]) -> tuple[Vertex, ...]:
    """Return all marked vertices in canonical order."""

    return tuple(sorted({vertex for edge in edges for vertex in edge}))


def enumerate_from_start(
    raw_start: Boundaries,
) -> tuple[frozenset[Boundaries], tuple[TerminalState, ...]]:
    """Enumerate every state and terminal leaf below one initial pattern."""

    start = canonicalize(raw_start)
    queue = deque((start,))
    seen = {start}
    leaves = []
    while queue:
        state = queue.popleft()
        children = transfer_children(state)
        if not children:
            leaves.append(terminal_state(state))
            continue
        for child in children:
            if child not in seen:
                seen.add(child)
                queue.append(child)
    return frozenset(seen), tuple(leaves)


def safe_decay(terminal: TerminalState) -> Fraction:
    """Return the accepted worst safe decay, capped at one full power."""

    audit = sigma_one_audit(terminal)
    if audit is None:
        return Fraction(1)
    return min(Fraction(1), audit.best_decay)


def full_image_audit() -> FullImageAudit:
    """Enumerate the complete terminal image from all initial collisions."""

    all_states: set[Boundaries] = set()
    terminal_by_boundaries: dict[Boundaries, TerminalState] = {}
    for start in INITIAL_BOUNDARIES:
        states, terminals = enumerate_from_start(start)
        all_states.update(states)
        for terminal in terminals:
            terminal_by_boundaries[terminal.boundaries] = terminal

    terminals = tuple(
        sorted(
            terminal_by_boundaries.values(),
            key=lambda terminal: (terminal.level, terminal.boundaries),
        )
    )
    sensitive = tuple(terminal for terminal in terminals if terminal.sensitive)
    counts = tuple(
        (level, sum(terminal.level == level for terminal in terminals))
        for level in range(4, 13)
    )
    sensitive_counts = tuple(
        (level, sum(terminal.level == level for terminal in sensitive))
        for level in range(4, 13)
    )
    minimum_decays = tuple(
        (
            level,
            min(
                safe_decay(terminal)
                for terminal in sensitive
                if terminal.level == level
            ),
        )
        for level in range(4, 13)
    )
    level_eight = tuple(terminal for terminal in sensitive if terminal.level == 8)
    level_seven_saturators = tuple(
        terminal
        for terminal in sensitive
        if terminal.level == 7 and safe_decay(terminal) == Fraction(1, 2)
    )
    level_six_saturators = tuple(
        terminal
        for terminal in sensitive
        if terminal.level == 6 and safe_decay(terminal) == Fraction(1, 2)
    )
    if tuple(terminal.boundaries for terminal in level_seven_saturators) != (
        LEVEL_SEVEN_BOUNDARIES,
    ):
        raise AssertionError(("level-seven saturators", level_seven_saturators))
    if tuple(terminal.boundaries for terminal in level_six_saturators) != (
        LEVEL_SIX_BOUNDARIES,
    ):
        raise AssertionError(("level-six saturators", level_six_saturators))

    return FullImageAudit(
        initial_potentials=tuple(
            branching_potential(canonicalize(start)) for start in INITIAL_BOUNDARIES
        ),
        reachable_states=len(all_states),
        terminal_types=len(terminals),
        sensitive_terminal_types=len(sensitive),
        terminal_types_by_level=counts,
        sensitive_types_by_level=sensitive_counts,
        minimum_safe_decay_by_level=minimum_decays,
        sensitive_level_eight_types=len(level_eight),
        level_eight_minimum_decay=min(map(safe_decay, level_eight)),
        level_seven_saturators=len(level_seven_saturators),
        level_six_saturators=len(level_six_saturators),
    )


def target_initial_states(
    target_boundaries: Boundaries,
) -> tuple[
    tuple[
        int,
        frozenset[Edge],
        dict[tuple[str, Vertex], LocalWeight],
    ],
    ...,
]:
    """Return every initial triple that can reach a fixed target."""

    target_edges = edges_from_boundaries(target_boundaries)
    target_level = len(vertices_from_edges(target_edges))
    boundary_edges = tuple(
        tuple(
            ((boundary_index, left), (boundary_index + 1, right))
            for left, right in edges
        )
        for boundary_index, edges in enumerate(target_boundaries)
    )
    result = []
    for chosen in product(*boundary_edges):
        edge_set = frozenset(chosen)
        potential = branching_potential(boundaries_from_edges(edge_set))
        if potential < target_level:
            continue
        weights: dict[tuple[str, Vertex], LocalWeight] = {}
        for left, right in chosen:
            weights[("A", left)] = LocalWeight("gamma")
            weights[("B", right)] = LocalWeight("gamma")
        result.append((potential, edge_set, weights))
    return tuple(result)


def target_transfer_options(
    target_edges: tuple[Edge, ...], edges: frozenset[Edge]
) -> tuple[tuple[str, Vertex, Edge], ...]:
    """Return every legal remaining edge toward a fixed terminal target."""

    outgoing, incoming = incidence(edges)
    remaining = set(target_edges) - set(edges)
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


def enumerate_target_histories(target_boundaries: Boundaries) -> CoefficientAudit:
    """Exhaust all legal initial triples and histories for one target."""

    target_edges = edges_from_boundaries(target_boundaries)
    target_edge_set = frozenset(target_edges)
    target_level = len(vertices_from_edges(target_edges))
    states = target_initial_states(target_boundaries)
    completed: list[
        tuple[
            int,
            tuple[TransferRecord, ...],
            tuple[tuple[str, Vertex, str, int], ...],
        ]
    ] = []
    contributing = Counter()

    def visit(
        initial_potential: int,
        edges: frozenset[Edge],
        weights: dict[tuple[str, Vertex], LocalWeight],
        records: tuple[TransferRecord, ...],
    ) -> None:
        if edges == target_edge_set:
            completed.append((initial_potential, records, weight_profile(weights)))
            return
        for direction, source, edge in target_transfer_options(target_edges, edges):
            new_edges, new_weights, record = apply_transfer(
                edges, weights, direction, source, edge
            )
            visit(
                initial_potential,
                new_edges,
                new_weights,
                records + (record,),
            )

    for potential, edges, weights in states:
        before = len(completed)
        visit(potential, edges, weights, ())
        if len(completed) > before:
            contributing[potential] += 1
    if not completed:
        raise AssertionError(("no histories", target_boundaries))

    initial_counts = Counter(potential for potential, _, _ in states)
    history_counts = Counter(potential for potential, _, _ in completed)
    profiles = Counter(profile for _, _, profile in completed)
    common_sites: set[tuple[str, Vertex]] | None = None
    weights_by_site: dict[tuple[str, Vertex], set[LocalWeight]] = {}
    for _, _, profile in completed:
        derivatives = [entry for entry in profile if entry[3] > 0]
        sites = {(side, vertex) for side, vertex, _, _ in derivatives}
        common_sites = sites if common_sites is None else common_sites & sites
        for side, vertex, kind, order in derivatives:
            weights_by_site.setdefault((side, vertex), set()).add(
                LocalWeight(kind, order)
            )
    if not common_sites:
        raise AssertionError(("no common derivative", target_boundaries))

    def transfer_counts(
        predicate: str,
    ) -> tuple[tuple[int, tuple[int, ...]], ...]:
        by_potential: dict[int, set[int]] = {}
        for potential, records, _ in completed:
            if predicate == "fresh":
                value = sum(record.creates_vertex for record in records)
            elif predicate == "existing":
                value = sum(not record.creates_vertex for record in records)
            elif predicate == "derivative":
                value = sum(record.differentiates_existing_weight for record in records)
            else:
                raise AssertionError(predicate)
            by_potential.setdefault(potential, set()).add(value)
        return tuple(
            (potential, tuple(sorted(values)))
            for potential, values in sorted(by_potential.items())
        )

    return CoefficientAudit(
        level=target_level,
        potential_initial_configurations=tuple(sorted(initial_counts.items())),
        contributing_initial_configurations=tuple(sorted(contributing.items())),
        histories_by_initial_potential=tuple(sorted(history_counts.items())),
        histories=len(completed),
        weight_profiles=len(profiles),
        fresh_transfers_by_initial_potential=transfer_counts("fresh"),
        existing_transfers_by_initial_potential=transfer_counts("existing"),
        derivative_events_by_initial_potential=transfer_counts("derivative"),
        common_derivative_sites=tuple(sorted(common_sites)),
        common_derivative_weights=tuple(
            (
                site,
                tuple(
                    sorted(
                        weights_by_site[site],
                        key=lambda weight: (weight.kind, weight.derivative),
                    )
                ),
            )
            for site in sorted(common_sites)
        ),
        time_exponents=tuple(
            len(boundary_edges) - 1 for boundary_edges in target_boundaries
        ),
    )


def dangerous_partitions(
    vertices: tuple[Vertex, ...], edges: tuple[Edge, ...]
) -> tuple[Partition, ...]:
    """Return all occupancy-two, cut-rank-two saturating partitions."""

    index = {vertex: position for position, vertex in enumerate(vertices)}
    indexed_edges = tuple((index[left], index[right]) for left, right in edges)
    result = []
    for partition in singleton_pair_partitions(tuple(range(len(vertices)))):
        if not any(len(block) == 2 for block in partition):
            continue
        if respecting_max_cut_rank(len(vertices), indexed_edges, partition) != 2:
            continue
        result.append(
            canonical_partition(
                [tuple(vertices[position] for position in block) for block in partition]
            )
        )
    return tuple(sorted(set(result)))


def low_level_repair_audit() -> LowLevelRepairAudit:
    """Score every retained level-seven and level-six repair branch."""

    level_seven_edges = edges_from_boundaries(LEVEL_SEVEN_BOUNDARIES)
    level_seven_vertices = vertices_from_edges(level_seven_edges)
    partitions_seven = dangerous_partitions(level_seven_vertices, level_seven_edges)
    if any(
        best_decay(level_seven_vertices, level_seven_edges, partition) != Fraction(1, 2)
        for partition in partitions_seven
    ):
        raise AssertionError("level-seven dangerous score")

    a0 = (0, 0)
    marked_b = ((1, 0), (1, 1))
    marked_c = ((2, 0), (2, 1))
    marked_d = ((3, 0), (3, 1))
    marked_neighbor_decays = [
        best_decay(
            level_seven_vertices,
            level_seven_edges + ((a0, neighbor),),
            partition,
        )
        for neighbor in marked_b
        for partition in partitions_seven
    ]

    new_b = (1, 2)
    one_new_partitions = extend_partitions(partitions_seven, new_b)
    existing_middle_decays = [
        best_decay(
            tuple(sorted(level_seven_vertices + (new_b,))),
            level_seven_edges + ((a0, new_b), (new_b, middle)),
            partition,
        )
        for middle in marked_c
        for partition in one_new_partitions
    ]

    new_c = (2, 2)
    two_new_partitions = extend_partitions(partitions_seven, new_b, new_c)
    existing_endpoint_decays = [
        best_decay(
            tuple(sorted(level_seven_vertices + (new_b, new_c))),
            level_seven_edges + ((a0, new_b), (new_b, new_c), (new_c, endpoint)),
            partition,
        )
        for endpoint in marked_d
        for partition in two_new_partitions
    ]

    new_d = (3, 2)
    three_new_partitions = extend_partitions(partitions_seven, new_b, new_c, new_d)
    fresh_edges = level_seven_edges + (
        (a0, new_b),
        (new_b, new_c),
        (new_c, new_d),
    )
    fresh_vertices = tuple(sorted(level_seven_vertices + (new_b, new_c, new_d)))
    fresh_endpoint_decays = [
        best_decay(fresh_vertices, fresh_edges, partition)
        for partition in three_new_partitions
    ]

    level_six_edges = edges_from_boundaries(LEVEL_SIX_BOUNDARIES)
    level_six_vertices = vertices_from_edges(level_six_edges)
    partitions_six = dangerous_partitions(level_six_vertices, level_six_edges)
    if any(
        best_decay(level_six_vertices, level_six_edges, partition) != Fraction(1, 2)
        for partition in partitions_six
    ):
        raise AssertionError("level-six dangerous score")
    b0 = (1, 0)
    marked_a = (0, 0)
    marked_outer_decays = [
        best_decay(
            level_six_vertices,
            level_six_edges + ((marked_a, b0),),
            partition,
        )
        for partition in partitions_six
    ]
    fresh_a = (0, 1)
    cancelled_partitions = extend_partitions(partitions_six, fresh_a)
    if sum(vertex[0] == 0 for vertex in level_six_vertices) != 1:
        raise AssertionError("level-six sensitivity")
    if sum(vertex[0] == 0 for vertex in level_six_vertices + (fresh_a,)) != 2:
        raise AssertionError("level-six fresh cancellation")

    return LowLevelRepairAudit(
        level_seven_dangerous_partitions=len(partitions_seven),
        level_seven_marked_neighbor_cases=len(marked_neighbor_decays),
        level_seven_existing_middle_cases=len(existing_middle_decays),
        level_seven_existing_endpoint_cases=len(existing_endpoint_decays),
        level_seven_fresh_endpoint_cases=len(fresh_endpoint_decays),
        level_seven_minimum_marked_neighbor_decay=min(marked_neighbor_decays),
        level_seven_minimum_existing_middle_decay=min(existing_middle_decays),
        level_seven_minimum_existing_endpoint_decay=min(existing_endpoint_decays),
        level_seven_minimum_fresh_endpoint_decay=min(fresh_endpoint_decays),
        level_six_dangerous_partitions=len(partitions_six),
        level_six_marked_outer_cases=len(marked_outer_decays),
        level_six_cancelled_fresh_outer_cases=len(cancelled_partitions),
        level_six_minimum_marked_outer_decay=min(marked_outer_decays),
        proved_global_exponent=Fraction(1, 12),
    )


def complete_audit() -> CompleteAudit:
    """Run the complete low-level certificate and global exponent handoff."""

    image = full_image_audit()
    level_seven_coefficients = enumerate_target_histories(LEVEL_SEVEN_BOUNDARIES)
    level_six_coefficients = enumerate_target_histories(LEVEL_SIX_BOUNDARIES)
    repair = low_level_repair_audit()

    # Protect the two previously accepted high-level centered handoffs before
    # claiming the global exponent.  The remaining high-level types already
    # have one full power from the complete best-of-two audit.
    if level_nine_repair_audit().proved_global_exponent != Fraction(1, 16):
        raise AssertionError("level-nine handoff")
    if level_ten_repair_audit().proved_global_exponent != Fraction(1, 18):
        raise AssertionError("level-ten handoff")

    post_repair_decays = dict(image.minimum_safe_decay_by_level)
    for level in (6, 7, 9, 10):
        post_repair_decays[level] = Fraction(1)
    exponent = min(decay / level for level, decay in post_repair_decays.items())
    if exponent != repair.proved_global_exponent:
        raise AssertionError(("global exponent", exponent, repair))
    return CompleteAudit(
        image=image,
        level_seven_coefficients=level_seven_coefficients,
        level_six_coefficients=level_six_coefficients,
        repair=repair,
    )


def main() -> None:
    audit = complete_audit()
    image = audit.image
    seven = audit.level_seven_coefficients
    six = audit.level_six_coefficients
    repair = audit.repair
    print(
        "low-level terminal centered repair: "
        f"reachable={image.reachable_states},"
        f"terminals={image.terminal_types},"
        f"sensitive={image.sensitive_terminal_types},"
        f"by_level={image.terminal_types_by_level},"
        f"safe_by_level={image.minimum_safe_decay_by_level},"
        f"level8_types={image.sensitive_level_eight_types},"
        f"level8_decay={image.level_eight_minimum_decay},"
        f"level7_histories={seven.histories},"
        f"level7_profiles={seven.weight_profiles},"
        f"level7_partitions={repair.level_seven_dangerous_partitions},"
        f"level7_counts=("
        f"{repair.level_seven_marked_neighbor_cases},"
        f"{repair.level_seven_existing_middle_cases},"
        f"{repair.level_seven_existing_endpoint_cases},"
        f"{repair.level_seven_fresh_endpoint_cases}),"
        f"level6_histories={six.histories},"
        f"level6_partitions={repair.level_six_dangerous_partitions},"
        f"level6_cancelled={repair.level_six_cancelled_fresh_outer_cases},"
        f"global_exponent={repair.proved_global_exponent}"
    )


if __name__ == "__main__":
    main()
