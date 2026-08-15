#!/usr/bin/env python3
"""Enumerate high-level terminal Stein graphs and audit two safe norm routes.

For levels nine through twelve the initial internal covariance endpoints must
be distinct, so the initial branching potential is twelve.  This module
enumerates the exact set-valued Stein dynamics up to layered relabeling.

Only placements with assigned suppression integer one can miss the desired
level-sensitive N^{-v/16} row.  For each such terminal graph we compute the
worst safe grouped-entry projective decay.  The final score is the better of
the global assigned and all-projective bounds; the two induction regimes are
never mixed.
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product


LAYERS = 4
INITIAL_POTENTIAL = 12

Vertex = tuple[int, int]
BoundaryEdge = tuple[int, int]
Boundaries = tuple[tuple[BoundaryEdge, ...], ...]
BlockPartition = tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class TerminalState:
    """One canonical terminal layered graph."""

    boundaries: Boundaries
    vertices: tuple[Vertex, ...]
    edges: tuple[tuple[Vertex, Vertex], ...]
    components: tuple[tuple[Vertex, ...], ...]
    level: int
    first_layer_vertices: int
    sensitive: bool
    total_surplus: int


@dataclass(frozen=True)
class ComponentNormAudit:
    """Cut-rank data for one connected graph component."""

    vertices: int
    edges: int
    surplus: int
    singleton_rank: int
    singleton_decay: Fraction
    paired_worst_rank: int
    paired_worst_decay: Fraction
    paired_witness: BlockPartition


@dataclass(frozen=True)
class SigmaOneTerminalAudit:
    """Best-of-two score for one dangerous assigned-sigma-one type."""

    level: int
    boundaries: Boundaries
    component_sizes: tuple[int, ...]
    component_surpluses: tuple[int, ...]
    strong_component: int | None
    projective_decay: Fraction
    assigned_decay: Fraction
    best_decay: Fraction
    required_decay: Fraction
    passes: bool
    paired_partition: BlockPartition | None


@dataclass(frozen=True)
class HighLevelTerminalAudit:
    """Complete high-level enumeration and best-of-two verdict."""

    reachable_states: int
    terminal_states: int
    high_level_terminals: int
    sensitive_high_level_terminals: int
    terminals_by_level: tuple[tuple[int, int], ...]
    sensitive_by_level: tuple[tuple[int, int], ...]
    sigma_one_types: int
    passing_sigma_one_types: int
    failing_sigma_one_types: int
    worst_best_decay: Fraction | None
    worst_required_decay: Fraction | None
    worst_safe_decay_by_level: tuple[tuple[int, Fraction], ...]
    previous_global_exponent: Fraction
    proved_global_exponent: Fraction
    failures: tuple[SigmaOneTerminalAudit, ...]


def layer_vertices(boundaries: Boundaries) -> tuple[tuple[int, ...], ...]:
    """Return the occupied coordinate labels in every layer."""

    occupied = [set() for _ in range(LAYERS)]
    for boundary, edges in enumerate(boundaries):
        for left, right in edges:
            occupied[boundary].add(left)
            occupied[boundary + 1].add(right)
    return tuple(tuple(sorted(values)) for values in occupied)


@lru_cache(maxsize=None)
def canonicalize(boundaries: Boundaries) -> Boundaries:
    """Canonicalize a layered graph under independent within-layer relabeling."""

    occupied = layer_vertices(boundaries)
    best: Boundaries | None = None
    permutation_lists = [tuple(permutations(values)) for values in occupied]
    for orders in product(*permutation_lists):
        maps = [{old: new for new, old in enumerate(order)} for order in orders]
        transformed = tuple(
            tuple(
                sorted(
                    (maps[boundary][left], maps[boundary + 1][right])
                    for left, right in edges
                )
            )
            for boundary, edges in enumerate(boundaries)
        )
        if best is None or transformed < best:
            best = transformed
    if best is None:
        raise AssertionError("a Stein state cannot be empty")
    return best


def incidence_sets(
    boundaries: Boundaries,
) -> tuple[tuple[frozenset[int], ...], tuple[frozenset[int], ...]]:
    """Return outgoing-right and incoming-left incidence sets by layer."""

    outgoing = [set() for _ in range(LAYERS)]
    incoming = [set() for _ in range(LAYERS)]
    for boundary, edges in enumerate(boundaries):
        for left, right in edges:
            outgoing[boundary].add(left)
            incoming[boundary + 1].add(right)
    return tuple(map(frozenset, outgoing)), tuple(map(frozenset, incoming))


def branching_potential(boundaries: Boundaries) -> int:
    """Return the exact four-layer Stein branching potential."""

    outgoing, incoming = incidence_sets(boundaries)
    vertex_count = sum(len(values) for values in layer_vertices(boundaries))
    mismatch = 0
    for layer in (1, 2):
        mismatch += (3 - layer) * len(incoming[layer] - outgoing[layer])
        mismatch += layer * len(outgoing[layer] - incoming[layer])
    return vertex_count + mismatch


def active_transfers(boundaries: Boundaries) -> tuple[tuple[int, int, str], ...]:
    """Return active ``(layer, source, direction)`` transfers."""

    outgoing, incoming = incidence_sets(boundaries)
    transfers = []
    for layer in (1, 2):
        transfers.extend(
            (layer, source, "left")
            for source in sorted(outgoing[layer] - incoming[layer])
        )
        transfers.extend(
            (layer, source, "right")
            for source in sorted(incoming[layer] - outgoing[layer])
        )
    return tuple(transfers)


def transfer_children(boundaries: Boundaries) -> tuple[Boundaries, ...]:
    """Apply every legal next transfer, including every existing or fresh neighbor."""

    old_potential = branching_potential(boundaries)
    occupied = layer_vertices(boundaries)
    children = set()
    for layer, source, direction in active_transfers(boundaries):
        neighbor_layer = layer - 1 if direction == "left" else layer + 1
        existing = occupied[neighbor_layer]
        choices = existing + (len(existing),)
        for neighbor in choices:
            creates_vertex = neighbor == len(existing)
            edge_lists = [list(edges) for edges in boundaries]
            if direction == "left":
                edge_lists[layer - 1].append((neighbor, source))
            else:
                edge_lists[layer].append((source, neighbor))
            raw = tuple(tuple(sorted(edges)) for edges in edge_lists)
            child = canonicalize(raw)
            new_potential = branching_potential(child)
            if creates_vertex and new_potential != old_potential:
                raise AssertionError(("fresh transfer potential", boundaries, child))
            if not creates_vertex and not new_potential < old_potential:
                raise AssertionError(("existing transfer potential", boundaries, child))
            children.add(child)
    return tuple(sorted(children))


def initial_state() -> Boundaries:
    """Return the unique high-level initial state with no matched internal endpoint."""

    # Boundary edges are layer 0--1, 1--2, and 2--3.  The two labels in each
    # internal layer are distinct; otherwise the initial potential is at most
    # eight and no terminal level nine or higher is possible.
    return canonicalize((((0, 0),), ((1, 0),), ((1, 0),)))


def graph_data(
    boundaries: Boundaries,
) -> tuple[tuple[Vertex, ...], tuple[tuple[Vertex, Vertex], ...]]:
    """Return canonical vertices and ordinary layered graph edges."""

    vertices = tuple(
        (layer, coordinate)
        for layer, values in enumerate(layer_vertices(boundaries))
        for coordinate in values
    )
    edges = tuple(
        ((boundary, left), (boundary + 1, right))
        for boundary, boundary_edges in enumerate(boundaries)
        for left, right in boundary_edges
    )
    return vertices, edges


def connected_components(
    vertices: tuple[Vertex, ...],
    edges: tuple[tuple[Vertex, Vertex], ...],
) -> tuple[tuple[Vertex, ...], ...]:
    """Return connected components in canonical vertex order."""

    adjacency: dict[Vertex, set[Vertex]] = defaultdict(set)
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    seen = set()
    components = []
    for start in vertices:
        if start in seen:
            continue
        queue = deque((start,))
        seen.add(start)
        component = set()
        while queue:
            vertex = queue.popleft()
            component.add(vertex)
            for neighbor in adjacency[vertex]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        components.append(tuple(sorted(component)))
    return tuple(components)


def terminal_state(boundaries: Boundaries) -> TerminalState:
    """Construct terminal graph metadata and verify the leaf invariant."""

    if active_transfers(boundaries):
        raise ValueError("state is not terminal")
    vertices, edges = graph_data(boundaries)
    components = connected_components(vertices, edges)
    total_surplus = len(edges) - len(vertices) + len(components)
    first_layer = len(layer_vertices(boundaries)[0])
    if branching_potential(boundaries) != len(vertices):
        raise AssertionError(("terminal potential", boundaries))
    return TerminalState(
        boundaries=boundaries,
        vertices=vertices,
        edges=edges,
        components=components,
        level=len(vertices),
        first_layer_vertices=first_layer,
        sensitive=first_layer % 2 == 1,
        total_surplus=total_surplus,
    )


def enumerate_terminal_states() -> tuple[int, tuple[TerminalState, ...]]:
    """Enumerate every terminal state reachable from the high-level initial leaf."""

    start = initial_state()
    if branching_potential(start) != INITIAL_POTENTIAL:
        raise AssertionError(("initial potential", start))
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
    return len(seen), tuple(
        sorted(leaves, key=lambda leaf: (leaf.level, leaf.boundaries))
    )


def gf2_rank(rows: tuple[tuple[int, ...], ...]) -> int:
    """Return exact binary matrix rank."""

    if not rows:
        return 0
    width = len(rows[0])
    values = [
        sum((entry & 1) << column for column, entry in enumerate(row)) for row in rows
    ]
    rank = 0
    for column in range(width):
        pivot = next(
            (
                index
                for index in range(rank, len(values))
                if values[index] >> column & 1
            ),
            None,
        )
        if pivot is None:
            continue
        values[rank], values[pivot] = values[pivot], values[rank]
        for index in range(len(values)):
            if index != rank and values[index] >> column & 1:
                values[index] ^= values[rank]
        rank += 1
    return rank


def indexed_component(
    terminal: TerminalState,
    component: tuple[Vertex, ...],
) -> tuple[int, tuple[tuple[int, int], ...]]:
    """Return component size and indexed simple edges."""

    index = {vertex: number for number, vertex in enumerate(component)}
    component_set = set(component)
    edges = tuple(
        sorted(
            (index[left], index[right])
            for left, right in terminal.edges
            if left in component_set
        )
    )
    return len(component), edges


def cut_rank(
    vertices: int,
    edges: tuple[tuple[int, int], ...],
    selected: frozenset[int],
) -> int:
    """Return binary adjacency rank across one component cut."""

    complement = tuple(vertex for vertex in range(vertices) if vertex not in selected)
    selected_tuple = tuple(sorted(selected))
    edge_parity: dict[tuple[int, int], int] = defaultdict(int)
    for left, right in edges:
        edge_parity[tuple(sorted((left, right)))] ^= 1
    rows = tuple(
        tuple(edge_parity[tuple(sorted((left, right)))] for right in complement)
        for left in selected_tuple
    )
    return gf2_rank(rows)


def maximum_cut_rank(vertices: int, edges: tuple[tuple[int, int], ...]) -> int:
    """Return maximum binary cut rank over all vertex bipartitions."""

    return max(
        cut_rank(vertices, edges, frozenset(selected))
        for mask in range(1, 1 << (vertices - 1))
        for selected in (
            tuple(vertex for vertex in range(vertices) if mask >> vertex & 1),
        )
    )


def singleton_pair_partitions(vertices: tuple[int, ...]) -> tuple[BlockPartition, ...]:
    """Enumerate partitions into singleton and pair blocks with at least one pair."""

    if not vertices:
        return ((),)
    first = vertices[0]
    rest = vertices[1:]
    result = []
    for suffix in singleton_pair_partitions(rest):
        result.append(((first,),) + suffix)
    for index, partner in enumerate(rest):
        remaining = rest[:index] + rest[index + 1 :]
        for suffix in singleton_pair_partitions(remaining):
            result.append(((first, partner),) + suffix)
    if len(vertices) == 1:
        return tuple(result)
    canonical = {
        tuple(
            sorted(
                (tuple(sorted(block)) for block in partition),
                key=lambda block: block[0],
            )
        )
        for partition in result
    }
    return tuple(sorted(canonical))


def respecting_max_cut_rank(
    vertices: int,
    edges: tuple[tuple[int, int], ...],
    partition: BlockPartition,
) -> int:
    """Return maximum cut rank among unions of physical-entry blocks."""

    ranks = []
    for mask in range(1, 1 << len(partition)):
        if mask == (1 << len(partition)) - 1:
            continue
        selected = frozenset(
            vertex
            for block_index, block in enumerate(partition)
            if mask >> block_index & 1
            for vertex in block
        )
        ranks.append(cut_rank(vertices, edges, selected))
    return max(ranks)


@lru_cache(maxsize=None)
def component_norm_audit(
    vertices: int,
    edges: tuple[tuple[int, int], ...],
) -> ComponentNormAudit:
    """Audit singleton and worst maximum-occupancy-two projective routes."""

    edge_count = len(edges)
    surplus = edge_count - vertices + 1
    singleton_rank = maximum_cut_rank(vertices, edges)
    singleton_decay = Fraction(surplus + singleton_rank - 1, 2)

    paired_candidates = []
    for partition in singleton_pair_partitions(tuple(range(vertices))):
        if not any(len(block) == 2 for block in partition):
            continue
        rank = respecting_max_cut_rank(vertices, edges, partition)
        paired_candidates.append((rank, partition))
    paired_rank, paired_witness = min(
        paired_candidates, key=lambda item: (item[0], item[1])
    )
    paired_decay = Fraction(surplus + paired_rank - 1, 2)
    return ComponentNormAudit(
        vertices=vertices,
        edges=edge_count,
        surplus=surplus,
        singleton_rank=singleton_rank,
        singleton_decay=singleton_decay,
        paired_worst_rank=paired_rank,
        paired_worst_decay=paired_decay,
        paired_witness=paired_witness,
    )


def sigma_one_audit(terminal: TerminalState) -> SigmaOneTerminalAudit | None:
    """Return the worst safe projective score if assigned sigma can equal one."""

    component_data = [
        indexed_component(terminal, component) for component in terminal.components
    ]
    norms = [
        component_norm_audit(vertices, edges) for vertices, edges in component_data
    ]
    total_surplus = sum(norm.surplus for norm in norms)
    if total_surplus > 1:
        return None

    strong_component: int | None = None
    paired_partition: BlockPartition | None = None
    if total_surplus == 1:
        projective_decay = sum((norm.singleton_decay for norm in norms), Fraction())
    else:
        candidates = []
        for index, norm in enumerate(norms):
            decay = norm.paired_worst_decay + sum(
                (
                    other.singleton_decay
                    for other_index, other in enumerate(norms)
                    if other_index != index
                ),
                Fraction(),
            )
            candidates.append((decay, index, norm.paired_witness))
        projective_decay, strong_component, paired_partition = min(candidates)

    assigned_decay = Fraction(1, 2)
    best_decay = max(assigned_decay, projective_decay)
    required = Fraction(terminal.level, 16)
    return SigmaOneTerminalAudit(
        level=terminal.level,
        boundaries=terminal.boundaries,
        component_sizes=tuple(norm.vertices for norm in norms),
        component_surpluses=tuple(norm.surplus for norm in norms),
        strong_component=strong_component,
        projective_decay=projective_decay,
        assigned_decay=assigned_decay,
        best_decay=best_decay,
        required_decay=required,
        passes=best_decay >= required,
        paired_partition=paired_partition,
    )


def high_level_terminal_audit() -> HighLevelTerminalAudit:
    """Run the complete level-nine-through-twelve best-of-two audit."""

    reachable, terminals = enumerate_terminal_states()
    high = tuple(terminal for terminal in terminals if 9 <= terminal.level <= 12)
    sensitive = tuple(terminal for terminal in high if terminal.sensitive)
    audits = tuple(
        audit
        for terminal in sensitive
        if (audit := sigma_one_audit(terminal)) is not None
    )
    failures = tuple(audit for audit in audits if not audit.passes)
    counts = tuple(
        (level, sum(terminal.level == level for terminal in high))
        for level in range(9, 13)
    )
    sensitive_counts = tuple(
        (level, sum(terminal.level == level for terminal in sensitive))
        for level in range(9, 13)
    )
    worst = (
        min(
            audits,
            key=lambda audit: (
                audit.best_decay - audit.required_decay,
                audit.boundaries,
            ),
        )
        if audits
        else None
    )
    worst_safe_by_level = []
    for level in range(9, 13):
        level_decays = []
        for terminal in sensitive:
            if terminal.level != level:
                continue
            terminal_audit = sigma_one_audit(terminal)
            if terminal_audit is None:
                # Assigned sigma is then at least two for every placement.
                level_decays.append(Fraction(1))
            else:
                # Placements with assigned sigma at least two retain one full
                # power, while the audit is the exact worst sigma-one score.
                level_decays.append(min(Fraction(1), terminal_audit.best_decay))
        worst_safe_by_level.append((level, min(level_decays)))

    exponent_candidates = [Fraction(1, 2 * level) for level in range(4, 9)]
    exponent_candidates.extend(decay / level for level, decay in worst_safe_by_level)
    proved_global_exponent = min(exponent_candidates)
    return HighLevelTerminalAudit(
        reachable_states=reachable,
        terminal_states=len(terminals),
        high_level_terminals=len(high),
        sensitive_high_level_terminals=len(sensitive),
        terminals_by_level=counts,
        sensitive_by_level=sensitive_counts,
        sigma_one_types=len(audits),
        passing_sigma_one_types=sum(audit.passes for audit in audits),
        failing_sigma_one_types=len(failures),
        worst_best_decay=worst.best_decay if worst else None,
        worst_required_decay=worst.required_decay if worst else None,
        worst_safe_decay_by_level=tuple(worst_safe_by_level),
        previous_global_exponent=Fraction(1, 24),
        proved_global_exponent=proved_global_exponent,
        failures=failures,
    )


def main() -> None:
    audit = high_level_terminal_audit()
    print(
        "high-level terminal best-of-two audit: "
        f"reachable={audit.reachable_states},"
        f"terminals={audit.terminal_states},"
        f"high={audit.high_level_terminals},"
        f"sensitive_high={audit.sensitive_high_level_terminals},"
        f"by_level={audit.terminals_by_level},"
        f"sensitive_by_level={audit.sensitive_by_level},"
        f"sigma_one={audit.sigma_one_types},"
        f"passing={audit.passing_sigma_one_types},"
        f"failing={audit.failing_sigma_one_types},"
        f"worst_best={audit.worst_best_decay},"
        f"worst_required={audit.worst_required_decay},"
        f"safe_by_level={audit.worst_safe_decay_by_level},"
        f"global_exponent={audit.proved_global_exponent}"
    )
    for failure in audit.failures[:10]:
        print(
            "joint saturator: "
            f"level={failure.level},"
            f"components={failure.component_sizes},"
            f"surpluses={failure.component_surpluses},"
            f"strong={failure.strong_component},"
            f"projective={failure.projective_decay},"
            f"required={failure.required_decay},"
            f"partition={failure.paired_partition},"
            f"boundaries={failure.boundaries}"
        )


if __name__ == "__main__":
    main()
