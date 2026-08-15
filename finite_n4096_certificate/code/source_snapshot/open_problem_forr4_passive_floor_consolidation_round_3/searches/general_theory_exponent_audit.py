#!/usr/bin/env python3
"""Exact exponent ledger for the repaired reverse-tree contraction.

The accepted dose audit charges a level-v diagram by D^v, not D^(2v).
This script records the graph-dependent half-powers of N that remain visible
before the proof replaces all of them by one N^(-1/2), checks the corrected
theorem ladder, and constructs a legal layered graph/placement family on
which both existing contraction branches retain only that one half-power.

The saturating family is a limitation of the current graph interface.  It
does not rule out a stronger theorem that uses interpolation coefficients,
physical-frame coupling, or cancellation between terminal diagrams.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


LEVELS = tuple(range(4, 13))


@dataclass(frozen=True)
class LayeredTreeWitness:
    vertices: int
    layers: tuple[int, ...]
    edges: tuple[tuple[int, int], ...]
    natural_cut_rank: int
    projective_sigma: int
    assigned_sigma: int
    singleton_placement: tuple[int, ...]
    assigned_placement: tuple[int, ...]


@dataclass(frozen=True)
class LevelLedger:
    level: int
    accepted_n_power: Fraction
    accepted_floor: Fraction
    n_one_sixteenth_power: Fraction
    n_one_sixteenth_extra: Fraction
    second_half_power: Fraction
    n_one_eighth_power: Fraction


def gf2_rank(rows: tuple[tuple[int, ...], ...]) -> int:
    """Return the exact binary rank of a zero-one matrix."""

    if not rows:
        return 0
    width = len(rows[0])
    values = [
        sum((entry & 1) << column for column, entry in enumerate(row))
        for row in rows
    ]
    rank = 0
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(values)) if values[index] >> column & 1),
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


def natural_cut_matrix(
    layers: tuple[int, ...], edges: tuple[tuple[int, int], ...]
) -> tuple[tuple[int, ...], ...]:
    """Return the layer-(1,3) versus layer-(2,4) adjacency over F_2."""

    left = tuple(vertex for vertex, layer in enumerate(layers) if layer in (1, 3))
    right = tuple(vertex for vertex, layer in enumerate(layers) if layer in (2, 4))
    left_index = {vertex: index for index, vertex in enumerate(left)}
    right_index = {vertex: index for index, vertex in enumerate(right)}
    matrix = [[0] * len(right) for _ in left]
    for first, second in edges:
        if first in left_index:
            row = left_index[first]
            column = right_index[second]
        else:
            row = left_index[second]
            column = right_index[first]
        matrix[row][column] ^= 1
    return tuple(tuple(row) for row in matrix)


def layered_tree_witness(vertices: int) -> LayeredTreeWitness:
    """Build a four-layer tree saturating the current suppression interface.

    Vertices 0--3 form one path through the four layers.  Every additional
    vertex is a layer-one leaf attached to the layer-two path vertex.  The
    One placement puts every vertex in its own physical amplitude entry and
    exercises the projective branch.  A second puts exactly two vertices in
    one entry and every other vertex in its own entry, exercising the
    assigned-fiber branch.
    """

    if vertices < 4:
        raise ValueError(("four-layer witness needs at least four vertices", vertices))
    layers = (1, 2, 3, 4) + (1,) * (vertices - 4)
    edges = ((0, 1), (1, 2), (2, 3)) + tuple(
        (vertex, 1) for vertex in range(4, vertices)
    )
    matrix = natural_cut_matrix(layers, edges)
    rank = gf2_rank(matrix)
    edge_count = len(edges)
    singleton_placement = tuple(range(vertices))
    assigned_placement = (0, 0) + tuple(range(1, vertices - 1))
    maximum_entry_occupancy = max(
        assigned_placement.count(entry) for entry in set(assigned_placement)
    )
    projective_sigma = edge_count + rank - vertices
    assigned_sigma = edge_count - vertices + maximum_entry_occupancy
    return LayeredTreeWitness(
        vertices=vertices,
        layers=layers,
        edges=edges,
        natural_cut_rank=rank,
        projective_sigma=projective_sigma,
        assigned_sigma=assigned_sigma,
        singleton_placement=singleton_placement,
        assigned_placement=assigned_placement,
    )


def accepted_floor_for_level(level: int) -> Fraction:
    """Floor exponent from D^level N^(-1/2)."""

    return Fraction(1, 2 * level)


def corrected_uniform_floor(
    n_power_per_vertex: Fraction,
    dose_power_per_vertex: int = 1,
) -> Fraction:
    """Floor exponent from D^(a v) N^(-b v)."""

    return n_power_per_vertex / dose_power_per_vertex


def level_ledger(level: int) -> LevelLedger:
    """Return accepted and target N-powers for one interpolation level."""

    accepted = Fraction(1, 2)
    one_sixteenth = max(accepted, Fraction(level, 16))
    second_half_power = Fraction(1) if level >= 9 else accepted
    return LevelLedger(
        level=level,
        accepted_n_power=accepted,
        accepted_floor=accepted / level,
        n_one_sixteenth_power=one_sixteenth,
        n_one_sixteenth_extra=one_sixteenth - accepted,
        second_half_power=second_half_power,
        n_one_eighth_power=Fraction(level, 8),
    )


def transcript_floor(
    n_powers: dict[int, Fraction],
) -> Fraction:
    """Return the floor exponent forced by all levels with D^v cost."""

    return min(n_powers[level] / level for level in LEVELS)


def accepted_transcript_floor() -> Fraction:
    return transcript_floor({level: Fraction(1, 2) for level in LEVELS})


def second_suppression_floor() -> Fraction:
    """Floor if levels nine through twelve gain a second N^-1/2."""

    return transcript_floor(
        {
            level: Fraction(1) if level >= 9 else Fraction(1, 2)
            for level in LEVELS
        }
    )


def one_eighth_floor() -> Fraction:
    return transcript_floor(
        {level: Fraction(level, 8) for level in LEVELS}
    )


def audit() -> tuple[tuple[LayeredTreeWitness, ...], tuple[LevelLedger, ...]]:
    witnesses = tuple(layered_tree_witness(level) for level in LEVELS)
    for witness in witnesses:
        vertices = witness.vertices
        if len(witness.edges) != vertices - 1:
            raise AssertionError(("witness is not a tree", witness))
        if set(witness.layers) != {1, 2, 3, 4}:
            raise AssertionError(("witness misses a layer", witness))
        if any(abs(witness.layers[u] - witness.layers[v]) != 1 for u, v in witness.edges):
            raise AssertionError(("nonadjacent layered edge", witness))
        if 4 * len(witness.edges) < 3 * vertices:
            raise AssertionError(("interpolation edge floor", witness))
        if witness.natural_cut_rank != 2:
            raise AssertionError(("natural-cut rank", witness))
        if witness.projective_sigma != 1 or witness.assigned_sigma != 1:
            raise AssertionError(("interface saturation", witness))
        if len(set(witness.singleton_placement)) != vertices:
            raise AssertionError(("projective placement is not singleton", witness))
        if max(
            witness.assigned_placement.count(entry)
            for entry in set(witness.assigned_placement)
        ) != 2:
            raise AssertionError(("assigned placement maximum", witness))

    ledgers = tuple(level_ledger(level) for level in LEVELS)
    if accepted_transcript_floor() != Fraction(1, 24):
        raise AssertionError("accepted exponent changed")
    if corrected_uniform_floor(Fraction(1, 8), 2) != Fraction(1, 16):
        raise AssertionError("historical duplicate-charge rung changed")
    if corrected_uniform_floor(Fraction(1, 8), 1) != Fraction(1, 8):
        raise AssertionError("corrected N^-v/8 rung changed")
    if second_suppression_floor() != Fraction(1, 16):
        raise AssertionError("second-suppression target changed")
    if one_eighth_floor() != Fraction(1, 8):
        raise AssertionError("one-eighth target changed")
    return witnesses, ledgers


def main() -> None:
    witnesses, ledgers = audit()
    print(
        "corrected exponent ladder: "
        f"accepted={accepted_transcript_floor()},"
        "stale_D2v_Nv8="
        f"{corrected_uniform_floor(Fraction(1, 8), 2)},"
        "current_Dv_Nv8="
        f"{corrected_uniform_floor(Fraction(1, 8), 1)},"
        f"second_suppression={second_suppression_floor()}"
    )
    for witness, ledger in zip(witnesses, ledgers, strict=True):
        print(
            "level audit: "
            f"v={witness.vertices},e={len(witness.edges)},"
            f"rank={witness.natural_cut_rank},"
            f"projective_sigma={witness.projective_sigma},"
            f"assigned_sigma={witness.assigned_sigma},"
            f"accepted_floor={ledger.accepted_floor},"
            f"N1_16_extra={ledger.n_one_sixteenth_extra},"
            f"N1_8_power={ledger.n_one_eighth_power}"
        )
    print("N=1024 second-half-power gain=1/32")


if __name__ == "__main__":
    main()
