#!/usr/bin/env python3
"""Sharpness audit for the limiting level-twelve three-path contraction.

The unique sensitive level-twelve terminal type is a forest of three
four-layer Hadamard paths.  A legal physical placement groups the first two
vertices of one path in a single amplitude entry and leaves every other
vertex singleton.  The accepted all-projective argument gives N^-1.

This module constructs unit grouped-entry test vectors whose contraction
with the *same-layer-distinct* graph tensor is

    (N - 1)^2 (N - 2)^3 / N^6
      = N^-1 (1 - N^-1)^2 (1 - 2 N^-1)^3.

Consequently the N^-1 dimension power is sharp for this graph-norm route.
Any asymptotic improvement beyond N^(1/12) must exploit additional frame
structure, cancel the positive coefficient family before taking norms, or
change the hard instance.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations
from math import isclose, sqrt

from high_level_terminal_best_of_two_audit import enumerate_terminal_states
from level_ten_forest_mean_zero_repair import (
    LocalWeight,
    TransferRecord,
    apply_transfer,
)
from low_level_terminal_centered_repair import (
    edges_from_boundaries,
    target_initial_states,
    target_transfer_options,
)
from terminal_interpolation_sigma_one_witness import TRANSFERS, replay_terminal_witness
from terminal_three_path_projective_repair import three_path_projective_certificate


def is_power_of_two(value: int) -> bool:
    """Return whether ``value`` is a positive power of two."""

    return value > 0 and value & (value - 1) == 0


def sylvester_entry(row: int, column: int, dimension: int) -> float:
    """Return one entry of the normalized Sylvester matrix."""

    if not is_power_of_two(dimension):
        raise ValueError(("Sylvester dimension", dimension))
    if not 0 <= row < dimension or not 0 <= column < dimension:
        raise ValueError(("Sylvester index", row, column, dimension))
    sign = -1.0 if (row & column).bit_count() % 2 else 1.0
    return sign / sqrt(dimension)


def masked_lower_witness_exact(dimension: int) -> Fraction:
    """Return the exact distinctness-masked lower-witness value."""

    if dimension < 3:
        raise ValueError(("three distinct coordinates required", dimension))
    return Fraction(
        (dimension - 1) ** 2 * (dimension - 2) ** 3,
        dimension**6,
    )


def lower_witness_ratio_to_n_inverse(dimension: int) -> Fraction:
    """Return the witness divided by ``N^-1``."""

    return dimension * masked_lower_witness_exact(dimension)


def direct_masked_contraction(dimension: int) -> float:
    """Contract the explicit unit vectors with every distinctness mask.

    The paths use vertices ``(a_j,b_j,c_j,d_j)``.  On path zero, ``a_0``
    and ``b_0`` form one grouped party with unit vector

        U(a,b) = H(a,b) H(b,r_0).

    Its other two parties are ``e_{r_0}`` and ``H e_{r_0}``.  On singleton
    path ``j``, use ``H e_{s_j}``, ``e_{s_j}``, ``e_{r_j}``, and a signed
    ``H e_{r_j}``.  Choose three distinct ``r`` values and two distinct
    ``s`` values.  The signs make every surviving summand nonnegative.

    The implementation sums the layerwise masked contraction directly.  It
    does not substitute the closed formula returned above.
    """

    if not is_power_of_two(dimension) or dimension < 4:
        raise ValueError(("power-of-two dimension at least four required", dimension))

    r0, r1, r2 = 0, 1, 2
    s1, s2 = 0, 1

    # The a-layer sum depends on b0 through U(a0,b0), so retain one value for
    # each allowed b0.  The a0/a1/a2 indices are required to be distinct.
    a_sums = {}
    for b0 in range(dimension):
        total = 0.0
        for a0 in range(dimension):
            paired = sylvester_entry(a0, b0, dimension) * sylvester_entry(
                b0, r0, dimension
            )
            paired *= sylvester_entry(a0, b0, dimension)
            for a1 in range(dimension):
                if a1 == a0:
                    continue
                first = sylvester_entry(a1, s1, dimension) ** 2
                for a2 in range(dimension):
                    if a2 == a0 or a2 == a1:
                        continue
                    second = sylvester_entry(a2, s2, dimension) ** 2
                    total += paired * first * second
        a_sums[b0] = total

    # b1=s1 and b2=s2 are fixed.  The b-layer mask leaves b0 outside those
    # two values.  The c coordinates r0,r1,r2 are distinct by construction.
    bc_sum = 0.0
    for b0 in range(dimension):
        if b0 in (s1, s2):
            continue
        bc_sum += (
            a_sums[b0]
            * sylvester_entry(b0, r0, dimension)
            * sylvester_entry(s1, r1, dimension)
            * sylvester_entry(s2, r2, dimension)
        )

    sign1 = 1.0 if sylvester_entry(s1, r1, dimension) > 0 else -1.0
    sign2 = 1.0 if sylvester_entry(s2, r2, dimension) > 0 else -1.0
    d_sum = 0.0
    for d0 in range(dimension):
        paired = sylvester_entry(r0, d0, dimension) ** 2
        for d1 in range(dimension):
            if d1 == d0:
                continue
            first = sign1 * sylvester_entry(r1, d1, dimension) ** 2
            for d2 in range(dimension):
                if d2 == d0 or d2 == d1:
                    continue
                second = sign2 * sylvester_entry(r2, d2, dimension) ** 2
                d_sum += paired * first * second
    return bc_sum * d_sum


def grouped_test_vector_norms(dimension: int) -> tuple[float, ...]:
    """Return the norms of all eleven grouped-entry test vectors."""

    if not is_power_of_two(dimension) or dimension < 4:
        raise ValueError(("power-of-two dimension at least four required", dimension))
    r0, r1, r2 = 0, 1, 2
    s1, s2 = 0, 1
    paired_squared = sum(
        (sylvester_entry(a, b, dimension) * sylvester_entry(b, r0, dimension)) ** 2
        for a in range(dimension)
        for b in range(dimension)
    )
    flat_squared = (
        sum(
            sylvester_entry(index, coordinate, dimension) ** 2
            for index in range(dimension)
            for coordinate in (r0, r1, r2, s1, s2)
        )
        / 5
    )
    # One paired vector, five basis vectors, and five Hadamard columns/rows.
    return (sqrt(paired_squared),) + (1.0,) * 5 + (sqrt(flat_squared),) * 5


def valid_positive_history_count() -> int:
    """Count legal orders for the previously displayed initial triple."""

    target = replay_terminal_witness()
    target_edges = set(target.edges)
    count = 0
    for order in permutations(TRANSFERS):
        try:
            witness = replay_terminal_witness(order)
        except ValueError:
            continue
        if witness.vertices != target.vertices or set(witness.edges) != target_edges:
            raise AssertionError(("valid order changes terminal support", order))
        if not witness.local_weight_strictly_positive:
            raise AssertionError(("valid order loses positivity", order))
        count += 1
    return count


def complete_positive_history_audit() -> tuple[int, int]:
    """Return all initial triples and histories for the canonical type.

    The earlier witness fixes one initial edge triple.  The coefficient of
    the canonical terminal support receives contributions from every legal
    initial triple contained in the three-path forest.  This recursion
    exhausts all of them and checks that every transfer is fresh and no
    existing positive local weight is differentiated.
    """

    _, terminals = enumerate_terminal_states()
    targets = [terminal for terminal in terminals if terminal.level == 12]
    if len(targets) != 1 or not targets[0].sensitive:
        raise AssertionError(("unique sensitive level-twelve type", targets))
    target_edges = edges_from_boundaries(targets[0].boundaries)
    target_edge_set = frozenset(target_edges)
    initial_states = target_initial_states(targets[0].boundaries)
    completed: list[tuple[TransferRecord, ...]] = []

    def visit(
        edges: frozenset[tuple[tuple[int, int], tuple[int, int]]],
        weights: dict[tuple[str, tuple[int, int]], LocalWeight],
        records: tuple[TransferRecord, ...],
    ) -> None:
        if edges == target_edge_set:
            completed.append(records)
            return
        for direction, source, edge in target_transfer_options(target_edges, edges):
            new_edges, new_weights, record = apply_transfer(
                edges, weights, direction, source, edge
            )
            visit(new_edges, new_weights, records + (record,))

    for potential, edges, weights in initial_states:
        if potential != 12:
            raise AssertionError(("level-twelve initial potential", potential))
        visit(edges, weights, ())
    if any(
        len(records) != 6
        or not all(record.creates_vertex for record in records)
        or any(record.differentiates_existing_weight for record in records)
        for records in completed
    ):
        raise AssertionError("nonpositive level-twelve coefficient history")
    return len(initial_states), len(completed)


@dataclass(frozen=True)
class LevelTwelveSharpnessAudit:
    """Exact obstruction certificate for the current norm architecture."""

    initial_configurations: int
    terminal_histories: int
    displayed_history_orders: int
    positive_coefficient_family: bool
    legal_assigned_sigma: int
    projective_upper_exponent: Fraction
    masked_lower_exponent: Fraction
    dimension: int
    masked_lower_value: Fraction
    ratio_to_n_inverse: Fraction


def sharpness_audit(dimension: int = 1024) -> LevelTwelveSharpnessAudit:
    """Return the level-twelve upper/lower exponent certificate."""

    witness = replay_terminal_witness()
    certificate = three_path_projective_certificate()
    initial_configurations, terminal_histories = complete_positive_history_audit()
    if certificate.combined_exponent != -1:
        raise AssertionError(("projective upper exponent", certificate))
    if not witness.local_weight_strictly_positive:
        raise AssertionError(("positive coefficient family", witness))
    lower = masked_lower_witness_exact(dimension)
    # The exact lower value is N^-1 times a factor tending to one, hence its
    # asymptotic dimension exponent is also -1.
    return LevelTwelveSharpnessAudit(
        initial_configurations=initial_configurations,
        terminal_histories=terminal_histories,
        displayed_history_orders=valid_positive_history_count(),
        positive_coefficient_family=witness.local_weight_strictly_positive,
        legal_assigned_sigma=witness.assigned_sigma,
        projective_upper_exponent=certificate.combined_exponent,
        masked_lower_exponent=Fraction(-1),
        dimension=dimension,
        masked_lower_value=lower,
        ratio_to_n_inverse=dimension * lower,
    )


def main() -> None:
    for dimension in (4, 8):
        direct = direct_masked_contraction(dimension)
        exact = float(masked_lower_witness_exact(dimension))
        if not isclose(direct, exact, rel_tol=0.0, abs_tol=2e-13):
            raise AssertionError(
                ("direct masked contraction", dimension, direct, exact)
            )
        if any(
            not isclose(norm, 1.0, abs_tol=2e-13)
            for norm in grouped_test_vector_norms(dimension)
        ):
            raise AssertionError(("grouped test-vector norms", dimension))
    audit = sharpness_audit()
    print(
        "level-twelve contraction sharpness: "
        f"initials={audit.initial_configurations},"
        f"histories={audit.terminal_histories},"
        f"displayed_orders={audit.displayed_history_orders},"
        f"positive={audit.positive_coefficient_family},"
        f"assigned_sigma={audit.legal_assigned_sigma},"
        f"upper_exponent={audit.projective_upper_exponent},"
        f"lower_exponent={audit.masked_lower_exponent},"
        f"N{audit.dimension}_lower={float(audit.masked_lower_value):.12g},"
        f"N{audit.dimension}_ratio={float(audit.ratio_to_n_inverse):.12g}"
    )


if __name__ == "__main__":
    main()
