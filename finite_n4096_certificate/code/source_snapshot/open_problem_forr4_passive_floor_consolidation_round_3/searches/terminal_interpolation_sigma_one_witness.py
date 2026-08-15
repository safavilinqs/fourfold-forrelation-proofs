#!/usr/bin/env python3
"""Exact terminal-interpolation witness for the Track B exponent gate.

The accepted k=4 Stein branching starts from one edge between each pair of
adjacent layers.  An active internal-layer vertex is transferred left or
right; the neighboring coordinate can be new or already marked.  This file
replays an explicit all-new transfer path that terminates at level twelve as
three disjoint four-layer paths.

The diagram is reflection-sensitive and occurs with nonzero positive local
weight.  A legal physical-entry placement pairs two vertices in one path and
leaves every other vertex singleton.  The exact assigned-fiber suppression
integer is then one.  Thus the sigma=1 obstruction occurs in the true
terminal interpolation image, not only at the relaxed graph interface.
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass


LAYERS = 4

Vertex = tuple[int, int]
Edge = tuple[Vertex, Vertex]


@dataclass(frozen=True)
class Transfer:
    """One accepted Stein transfer from an active internal vertex."""

    direction: str
    layer: int
    source: int
    neighbor: int
    creates_vertex: bool


@dataclass(frozen=True)
class TerminalInterpolationWitness:
    """Exact combinatorial and suppression data for one terminal path."""

    vertices: tuple[Vertex, ...]
    edges: tuple[Edge, ...]
    components: tuple[tuple[Vertex, ...], ...]
    initial_potential: int
    terminal_potential: int
    new_vertex_transfers: int
    existing_vertex_transfers: int
    first_layer_vertices: int
    reflection_sensitive: bool
    projective_sigma: int
    assigned_sigma: int
    maximum_boundary_degree: int
    paired_entry: tuple[Vertex, Vertex]
    local_weight_strictly_positive: bool


INITIAL_EDGES: tuple[Edge, ...] = (
    ((0, 0), (1, 0)),
    ((1, 3), (2, 2)),
    ((2, 1), (3, 1)),
)

# The two initial internal endpoints differ in both internal layers.  Each
# transfer below chooses a fresh neighbor.  Transfers two and four form one
# path, transfers one and five form another, and transfers three and six form
# the third.
TRANSFERS: tuple[Transfer, ...] = (
    Transfer("left", 1, 3, 1, True),
    Transfer("right", 1, 0, 0, True),
    Transfer("left", 2, 1, 2, True),
    Transfer("right", 2, 0, 0, True),
    Transfer("right", 2, 2, 2, True),
    Transfer("left", 1, 2, 3, True),
)


def marked_vertices(
    left_incidence: list[set[int]],
    right_incidence: list[set[int]],
) -> set[Vertex]:
    """Return all currently marked layer/coordinate pairs."""

    return {
        (layer, coordinate)
        for layer in range(LAYERS)
        for coordinate in left_incidence[layer] | right_incidence[layer]
    }


def branching_potential(
    left_incidence: list[set[int]],
    right_incidence: list[set[int]],
) -> int:
    """Return the exact k=4 Bansal--Sinha branching potential."""

    union_size = sum(
        len(left_incidence[layer] | right_incidence[layer]) for layer in range(LAYERS)
    )
    mismatch = 0
    for layer in (1, 2):
        mismatch += (3 - layer) * len(right_incidence[layer] - left_incidence[layer])
        mismatch += layer * len(left_incidence[layer] - right_incidence[layer])
    return union_size + mismatch


def connected_components(
    vertices: set[Vertex],
    edges: tuple[Edge, ...] | list[Edge],
) -> tuple[tuple[Vertex, ...], ...]:
    """Return the graph components in canonical order."""

    adjacency: dict[Vertex, set[Vertex]] = defaultdict(set)
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    seen: set[Vertex] = set()
    result = []
    for start in sorted(vertices):
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
        result.append(tuple(sorted(component)))
    return tuple(result)


def gf2_rank(rows: tuple[tuple[int, ...], ...]) -> int:
    """Return the exact binary rank of a zero-one matrix."""

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


def natural_cut_rank(component: tuple[Vertex, ...], edges: tuple[Edge, ...]) -> int:
    """Return the layer-(1,3) versus layer-(2,4) rank of one component."""

    component_set = set(component)
    left = tuple(vertex for vertex in component if vertex[0] in (0, 2))
    right = tuple(vertex for vertex in component if vertex[0] in (1, 3))
    left_index = {vertex: index for index, vertex in enumerate(left)}
    right_index = {vertex: index for index, vertex in enumerate(right)}
    matrix = [[0] * len(right) for _ in left]
    for first, second in edges:
        if first not in component_set:
            continue
        if first in left_index:
            row = left_index[first]
            column = right_index[second]
        else:
            row = left_index[second]
            column = right_index[first]
        matrix[row][column] ^= 1
    return gf2_rank(tuple(tuple(row) for row in matrix))


def component_edge_count(component: tuple[Vertex, ...], edges: tuple[Edge, ...]) -> int:
    """Count edges internal to one connected component."""

    component_set = set(component)
    return sum(
        first in component_set and second in component_set for first, second in edges
    )


def physical_entry_placement(
    vertices: tuple[Vertex, ...],
) -> tuple[dict[Vertex, int], tuple[Vertex, Vertex]]:
    """Pair two marks in the first path and make all others singleton."""

    paired = ((0, 0), (1, 0))
    placement = {paired[0]: 0, paired[1]: 0}
    next_entry = 1
    for vertex in vertices:
        if vertex in placement:
            continue
        placement[vertex] = next_entry
        next_entry += 1
    return placement, paired


def suppression_parameters(
    components: tuple[tuple[Vertex, ...], ...],
    edges: tuple[Edge, ...],
    placement: dict[Vertex, int],
) -> tuple[int, int]:
    """Return the exact projective and assigned-fiber sigma integers."""

    projective = 0
    assigned = 0
    for component in components:
        edge_count = component_edge_count(component, edges)
        surplus = edge_count - len(component) + 1
        rank = natural_cut_rank(component, edges)
        occupancy: dict[int, int] = defaultdict(int)
        for vertex in component:
            occupancy[placement[vertex]] += 1
        maximum_occupancy = max(occupancy.values())
        projective += surplus + rank - 1
        assigned += surplus + maximum_occupancy - 1
    return projective, assigned


def replay_terminal_witness(
    transfers: tuple[Transfer, ...] = TRANSFERS,
) -> TerminalInterpolationWitness:
    """Replay the explicit path and return its exact terminal data."""

    left_incidence = [set() for _ in range(LAYERS)]
    right_incidence = [set() for _ in range(LAYERS)]
    edges: list[Edge] = []
    for first, second in INITIAL_EDGES:
        if second[0] != first[0] + 1:
            raise ValueError(("initial edge is not adjacent", first, second))
        left_incidence[first[0]].add(first[1])
        right_incidence[second[0]].add(second[1])
        edges.append((first, second))

    initial_potential = branching_potential(left_incidence, right_incidence)
    new_transfers = 0
    existing_transfers = 0
    for transfer in transfers:
        if transfer.direction not in {"left", "right"}:
            raise ValueError(("transfer direction", transfer))
        layer = transfer.layer
        if layer not in (1, 2):
            raise ValueError(("transfer layer", transfer))
        source = transfer.source
        before = marked_vertices(left_incidence, right_incidence)
        before_potential = branching_potential(left_incidence, right_incidence)
        if transfer.direction == "right":
            if source not in right_incidence[layer] - left_incidence[layer]:
                raise ValueError(("inactive right transfer", transfer))
            neighbor = (layer + 1, transfer.neighbor)
            left_incidence[layer].add(source)
            right_incidence[layer + 1].add(transfer.neighbor)
            edge = ((layer, source), neighbor)
        else:
            if source not in left_incidence[layer] - right_incidence[layer]:
                raise ValueError(("inactive left transfer", transfer))
            neighbor = (layer - 1, transfer.neighbor)
            right_incidence[layer].add(source)
            left_incidence[layer - 1].add(transfer.neighbor)
            edge = (neighbor, (layer, source))
        created = neighbor not in before
        if created != transfer.creates_vertex:
            raise ValueError(("transfer growth flag", transfer, created))
        edges.append(edge)
        after_potential = branching_potential(left_incidence, right_incidence)
        if created:
            if after_potential != before_potential:
                raise AssertionError(("new-vertex potential", transfer))
            new_transfers += 1
        else:
            if not after_potential < before_potential:
                raise AssertionError(("existing-vertex potential", transfer))
            existing_transfers += 1

    for layer in (1, 2):
        if left_incidence[layer] != right_incidence[layer]:
            raise AssertionError(("nonterminal internal layer", layer))
    vertices = tuple(sorted(marked_vertices(left_incidence, right_incidence)))
    edge_tuple = tuple(edges)
    components = connected_components(set(vertices), edge_tuple)
    placement, paired = physical_entry_placement(vertices)
    projective_sigma, assigned_sigma = suppression_parameters(
        components,
        edge_tuple,
        placement,
    )
    boundary_degrees: dict[Vertex, int] = defaultdict(int)
    for first, second in edge_tuple:
        if first[0] == 0 and second[0] == 1:
            boundary_degrees[second] += 1
        if first[0] == 2 and second[0] == 3:
            boundary_degrees[first] += 1
    first_layer_vertices = sum(layer == 0 for layer, _ in vertices)

    # Every transfer creates a vertex.  Therefore no derivative lands on an
    # existing local weight.  The only scalar factors are positive Stein
    # kernels and psi' factors.  For the capped odd increasing psi,
    # S_psi(x)>0 and psi'(x)>0 for every finite x, so this path has a
    # pointwise strictly positive local weight and cannot disappear by a
    # scalar branching-sign cancellation.
    local_weight_positive = existing_transfers == 0

    return TerminalInterpolationWitness(
        vertices=vertices,
        edges=edge_tuple,
        components=components,
        initial_potential=initial_potential,
        terminal_potential=branching_potential(left_incidence, right_incidence),
        new_vertex_transfers=new_transfers,
        existing_vertex_transfers=existing_transfers,
        first_layer_vertices=first_layer_vertices,
        reflection_sensitive=first_layer_vertices % 2 == 1,
        projective_sigma=projective_sigma,
        assigned_sigma=assigned_sigma,
        maximum_boundary_degree=max(boundary_degrees.values()),
        paired_entry=paired,
        local_weight_strictly_positive=local_weight_positive,
    )


def relaxed_star_boundary_degree(vertices: int) -> int:
    """Return the layer-one degree of the old relaxed star witness."""

    if vertices < 4:
        raise ValueError(("four-layer star needs four vertices", vertices))
    return vertices - 3


def terminal_boundary_degree_cap() -> int:
    """Return the exact outer-boundary degree cap in the branching image.

    A layer-two vertex can receive one initial edge from layer one.  If it is
    unmatched toward layer one, it can be the source of exactly one leftward
    transfer, after which that mismatch is resolved.  No transfer is sourced
    in the outer layer.  Hence it has at most two layer-one neighbors.  The
    layer-three/layer-four boundary is symmetric.
    """

    return 2


def main() -> None:
    witness = replay_terminal_witness()
    print(
        "terminal interpolation sigma-one witness: "
        f"v={len(witness.vertices)},"
        f"e={len(witness.edges)},"
        f"components={len(witness.components)},"
        f"component_sizes={tuple(map(len, witness.components))},"
        f"initial_potential={witness.initial_potential},"
        f"terminal_potential={witness.terminal_potential},"
        f"new_transfers={witness.new_vertex_transfers},"
        f"existing_transfers={witness.existing_vertex_transfers},"
        f"first_layer={witness.first_layer_vertices},"
        f"sensitive={witness.reflection_sensitive},"
        f"projective_sigma={witness.projective_sigma},"
        f"assigned_sigma={witness.assigned_sigma},"
        f"max_boundary_degree={witness.maximum_boundary_degree},"
        f"terminal_boundary_cap={terminal_boundary_degree_cap()},"
        f"old_star_boundary_degree={relaxed_star_boundary_degree(12)},"
        f"positive_local_weight={witness.local_weight_strictly_positive}"
    )


if __name__ == "__main__":
    main()
