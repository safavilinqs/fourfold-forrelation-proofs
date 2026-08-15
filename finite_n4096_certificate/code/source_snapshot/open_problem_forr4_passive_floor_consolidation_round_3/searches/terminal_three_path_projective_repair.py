#!/usr/bin/env python3
"""All-projective repair for the true terminal three-path witness.

The terminal sigma-one witness consists of three disjoint four-layer paths.
One path has maximum physical-entry occupancy two, while the other two have
occupancy one.  The global assigned dichotomy retains only N^{-1/2}.

For this placement it is safe to keep the *entire* reverse frame skeleton in
the grouped-entry projective norm.  A four-layer path whose vertices occupy
at least two physical entries has grouped injective norm at most one.  If all
four vertices occupy distinct entries, its natural cut has binary rank two
and the norm is at most N^{-1/2}.  Vertical multiplicativity therefore gives
N^{-1} for the three-path forest.  No Hilbert auxiliary is reinterpreted as
a projective tensor, so the RT-003 counterexample is irrelevant.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product


PATH_VERTICES = (0, 1, 2, 3)
PATH_EDGES = ((0, 1), (1, 2), (2, 3))

Partition = tuple[tuple[int, ...], ...]


def canonical_partition(labels: tuple[int, ...]) -> Partition:
    """Return the set partition encoded by arbitrary integer labels."""

    relabel: dict[int, int] = {}
    blocks: list[list[int]] = []
    for vertex, label in enumerate(labels):
        if label not in relabel:
            relabel[label] = len(blocks)
            blocks.append([])
        blocks[relabel[label]].append(vertex)
    return tuple(tuple(block) for block in blocks)


def path_partitions() -> tuple[Partition, ...]:
    """Enumerate the fifteen physical-entry partitions of a four-vertex path."""

    return tuple(
        sorted(
            {
                canonical_partition(labels)
                for labels in product(
                    range(len(PATH_VERTICES)), repeat=len(PATH_VERTICES)
                )
            },
            key=lambda partition: (len(partition), partition),
        )
    )


def gf2_rank(rows: tuple[tuple[int, ...], ...]) -> int:
    """Return the exact rank of a zero-one matrix over GF(2)."""

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


def crossing_rank(vertices: frozenset[int]) -> int:
    """Return the path adjacency rank across ``vertices`` and its complement."""

    complement = tuple(vertex for vertex in PATH_VERTICES if vertex not in vertices)
    selected = tuple(sorted(vertices))
    matrix = []
    edge_set = {tuple(sorted(edge)) for edge in PATH_EDGES}
    for left in selected:
        matrix.append(
            tuple(int(tuple(sorted((left, right))) in edge_set) for right in complement)
        )
    return gf2_rank(tuple(matrix))


def flattening_exponent(vertices: frozenset[int]) -> Fraction:
    """Return ``a`` in the exact path flattening norm ``N**a``."""

    if not vertices or len(vertices) == len(PATH_VERTICES):
        raise ValueError(("flattening cut must be nontrivial", vertices))
    vertices_count = len(PATH_VERTICES)
    edges_count = len(PATH_EDGES)
    return Fraction(vertices_count - edges_count - crossing_rank(vertices), 2)


def grouped_injective_upper_exponent(partition: Partition) -> Fraction:
    """Best safe exponent obtained from cuts that respect physical entries.

    If all vertices occupy one entry, the grouped tensor is one Hilbert
    vector and its norm is the Frobenius norm N^{(v-e)/2}=N^{1/2}.
    Otherwise every union of entry blocks defines a legal flattening.  The
    injective norm is at most the smallest such operator norm.
    """

    if len(partition) == 1:
        return Fraction(len(PATH_VERTICES) - len(PATH_EDGES), 2)
    candidates = []
    for mask in range(1, 1 << len(partition)):
        if mask == (1 << len(partition)) - 1:
            continue
        selected = frozenset(
            vertex
            for block_index, block in enumerate(partition)
            if mask >> block_index & 1
            for vertex in block
        )
        candidates.append(flattening_exponent(selected))
    return min(candidates)


def maximum_occupancy(partition: Partition) -> int:
    """Return the largest number of path vertices in one physical entry."""

    return max(map(len, partition))


@dataclass(frozen=True)
class ThreePathProjectiveCertificate:
    """Exact exponent certificate for the terminal three-path forest."""

    partition_count: int
    strong_worst_exponent: Fraction
    weak_exponent: Fraction
    combined_exponent: Fraction
    accepted_exponent: Fraction
    level_twelve_target_exponent: Fraction
    distinctness_pairs: int

    @property
    def extra_gain_beyond_accepted(self) -> Fraction:
        return self.combined_exponent - self.accepted_exponent

    @property
    def target_slack(self) -> Fraction:
        return self.combined_exponent - self.level_twelve_target_exponent


def three_path_projective_certificate() -> ThreePathProjectiveCertificate:
    """Return the all-projective exponent certificate for occupancies (2,1,1)."""

    partitions = path_partitions()
    strong = [
        grouped_injective_upper_exponent(partition)
        for partition in partitions
        if maximum_occupancy(partition) == 2
    ]
    singleton = ((0,), (1,), (2,), (3,))
    weak = grouped_injective_upper_exponent(singleton)
    strong_worst = max(strong)
    combined = strong_worst + 2 * weak
    return ThreePathProjectiveCertificate(
        partition_count=len(partitions),
        strong_worst_exponent=strong_worst,
        weak_exponent=weak,
        combined_exponent=combined,
        accepted_exponent=Fraction(-1, 2),
        level_twelve_target_exponent=Fraction(-3, 4),
        distinctness_pairs=4 * 3,
    )


def main() -> None:
    certificate = three_path_projective_certificate()
    dimension = 1024
    print(
        "terminal three-path projective repair: "
        f"partitions={certificate.partition_count},"
        f"strong_worst={certificate.strong_worst_exponent},"
        f"weak={certificate.weak_exponent},"
        f"combined={certificate.combined_exponent},"
        f"accepted={certificate.accepted_exponent},"
        f"target={certificate.level_twelve_target_exponent},"
        f"extra_gain={certificate.extra_gain_beyond_accepted},"
        f"target_slack={certificate.target_slack},"
        f"mask_terms_at_most={2**certificate.distinctness_pairs},"
        f"N1024_bound={dimension ** float(certificate.combined_exponent):.12g},"
        f"N1024_target={dimension ** float(certificate.level_twelve_target_exponent):.12g}"
    )


if __name__ == "__main__":
    main()
