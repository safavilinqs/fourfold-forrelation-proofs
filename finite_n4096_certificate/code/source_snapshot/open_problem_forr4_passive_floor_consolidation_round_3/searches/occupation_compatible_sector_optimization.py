#!/usr/bin/env python3
"""Occupation-compatible one-batch sector certificate at N=1024.

For a profile ``a`` and occurrence split ``s``, a physical matrix entry with
row occupation ``n`` and column occupation ``m`` necessarily satisfies

    n - s = m - (a-s),

because both sides contain the same unmarked intersection support.  The
current moment relaxation replaces this exact pairing by a Cauchy product of
two totals over all occupations.  Here we retain the pairing and optimize the
resulting concave sum over the 210 dose-six occupation states.

The accepted partial ledger includes all profiles through degree eight, the
proved high endpoint profiles, the four leading triple-cubic profiles, and
both separated cubic--quintic reversal pairs.  Other degree-ten/twelve
profiles and the adaptive lift remain open.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from math import comb, prod, sqrt

import numpy as np

from attenuation_promise_concentration import promise_concentration
from degree_six_joint_occupation_optimization import (
    DEGREE_SIX_PROFILES,
    degree_six_coefficient,
)
from double_endpoint_occupation_optimization import (
    PROFILE as DOUBLE_ENDPOINT_PROFILE,
    coefficient as double_endpoint_coefficient,
    occupation_states,
)


ORDER = 32
BETA = 781 / 1000
MINIMAL_PROFILE = (1, 1, 1, 1)
DEGREE_EIGHT_PROFILES = tuple(
    profile
    for profile in product((1, 3, 5), repeat=4)
    if sum(profile) == 8
)
HIGH_DEGREE_PROFILES = tuple(
    profile
    for profile in product((1, 3, 5, 7, 9), repeat=4)
    if sum(profile) in (10, 12)
)
TRIPLE_CUBIC_PROFILES = (
    (3, 3, 3, 1),
    (1, 3, 3, 3),
    (3, 1, 3, 3),
    (3, 3, 1, 3),
)
SEPARATED_QUINTIC_CUBIC_PROFILES = (
    (5, 1, 3, 1),
    (1, 3, 1, 5),
)
SEPARATED_CUBIC_QUINTIC_PROFILES = (
    (3, 1, 5, 1),
    (1, 5, 1, 3),
)
KNOWN_HIGH_DEGREE_PROFILES = tuple(
    profile
    for profile in HIGH_DEGREE_PROFILES
    if (
        profile in TRIPLE_CUBIC_PROFILES
        or profile in SEPARATED_QUINTIC_CUBIC_PROFILES
        or profile in SEPARATED_CUBIC_QUINTIC_PROFILES
        or (
            len(
                tuple(
                    block
                    for block, degree in enumerate(profile)
                    if degree > 1
                )
            )
            == 1
            and next(
                block
                for block, degree in enumerate(profile)
                if degree > 1
            )
            in (0, 3)
        )
    )
)
TV_THRESHOLD = 1 / 3


@dataclass(frozen=True)
class CompatibleOccupationResult:
    objective: float
    supporting_upper: float
    profile_contributions: tuple[tuple[tuple[int, ...], float], ...]
    occupation_weights: tuple[tuple[tuple[int, ...], float], ...]
    support: tuple[tuple[float, tuple[int, ...]], ...]
    maximum_gradient_state: tuple[int, ...]
    leading_terms: tuple[
        tuple[
            float,
            tuple[int, ...],
            tuple[int, ...],
            tuple[int, ...],
            tuple[int, ...],
            float,
        ],
        ...,
    ]


def profile_splits(profile: tuple[int, ...]) -> list[tuple[int, ...]]:
    return list(product(*(range(degree + 1) for degree in profile)))


def path_cut_coefficient(split: tuple[int, ...]) -> float:
    mask = frozenset(block for block, selected in enumerate(split) if selected)
    complement = frozenset(range(4)) - mask
    canonical = min(
        mask, complement, key=lambda value: (len(value), sorted(value))
    )
    if not canonical:
        exponent = 3
    elif len(canonical) == 1:
        exponent = 2
    else:
        exponent = 2 if canonical in (
            frozenset({0, 1}),
            frozenset({2, 3}),
        ) else 1
    return ORDER ** (-exponent)


def coefficient(
    profile: tuple[int, ...],
    split: tuple[int, ...],
    include_record_three: bool,
    high_degree_coefficient: float | None,
) -> float:
    if profile == MINIMAL_PROFILE:
        return path_cut_coefficient(split)
    if profile == DOUBLE_ENDPOINT_PROFILE:
        return double_endpoint_coefficient(split)
    if sum(profile) == 8:
        return degree_eight_coefficient(profile, split, include_record_three)
    if sum(profile) > 8:
        if profile in TRIPLE_CUBIC_PROFILES:
            return triple_cubic_coefficient(profile, split)
        if profile in SEPARATED_QUINTIC_CUBIC_PROFILES:
            return separated_quintic_cubic_coefficient(profile, split)
        if profile in SEPARATED_CUBIC_QUINTIC_PROFILES:
            return separated_cubic_quintic_coefficient(profile, split)
        decorated = tuple(
            block for block, degree in enumerate(profile) if degree > 1
        )
        if len(decorated) == 1 and decorated[0] in (0, 3):
            smaller_side = min(sum(split), sum(profile) - sum(split))
            return min(1 / ORDER, ORDER ** (smaller_side - 3))
        if high_degree_coefficient is None:
            raise ValueError(("missing high-degree coefficient", profile))
        return high_degree_coefficient
    if profile in ((1, 3, 1, 1), (1, 1, 3, 1)):
        decorated = 1 if profile[1] == 3 else 2
        families = tuple(
            l_shape_incidence if block == decorated else singleton_incidence
            for block in range(4)
        )
        incidence = incidence_coefficient(
            profile,
            split,
            families,
            1 / (ORDER**3 * (ORDER - 1) ** 2),
        )
        return min(degree_six_coefficient(profile, split), incidence)
    return degree_six_coefficient(profile, split)


def adjacent_record_one_coefficient(
    endpoint_selected: int, middle_selected: int
) -> float:
    """Uniform-incidence bound (2.4) from the Round 2 adjacent note."""

    q = ORDER
    endpoint_degrees = (
        q * comb(q, 3) + q * q * (q - 1) * comb(q, 2),
        comb(q - 1, 2) + (q - 1) * comb(q, 2) + q * (q - 1) ** 2,
        q * q - 2,
        1,
    )
    middle_degrees = (
        q * q * (q - 1) ** 2,
        3 * (q - 1) ** 2,
        2 * (q - 1),
        1,
    )
    maximum_squared_entry = (q + 2) ** 2 / (
        q * q * (q - 1) ** 2 * (q - 2) ** 2
    )
    row = sqrt(
        maximum_squared_entry
        * endpoint_degrees[endpoint_selected]
        * middle_degrees[middle_selected]
    ) / (q - 1)
    column = sqrt(
        maximum_squared_entry
        * endpoint_degrees[3 - endpoint_selected]
        * middle_degrees[3 - middle_selected]
    ) / (q * (q - 1))
    return min(row, column)


def singleton_incidence(selected: int) -> int:
    return ORDER * ORDER if selected == 0 else 1


def endpoint_record_one_incidence(selected: int) -> int:
    """Cubic endpoint supports with one odd hidden label."""

    q = ORDER
    values = (
        q * comb(q, 3) + q * q * (q - 1) * comb(q, 2),
        comb(q - 1, 2)
        + (q - 1) * comb(q, 2)
        + q * (q - 1) ** 2,
        q * q - 2,
        1,
    )
    return values[selected]


def l_shape_incidence(selected: int) -> int:
    q = ORDER
    values = (q * q * (q - 1) ** 2, 3 * (q - 1) ** 2, 2 * (q - 1), 1)
    return values[selected]


def integer_partitions(total: int, maximum: int | None = None):
    """Yield decreasing positive partitions of a small integer."""

    if total == 0:
        yield ()
        return
    if maximum is None or maximum > total:
        maximum = total
    for first in range(maximum, 0, -1):
        for tail in integer_partitions(total - first, first):
            yield (first,) + tail


def one_record_extension_count(
    order: int, fixed_column_counts: tuple[int, ...], degree: int
) -> int:
    """Count degree-``degree`` supports containing a fixed column pattern."""

    fixed_size = sum(fixed_column_counts)
    columns = fixed_column_counts + (0,) * (
        order - len(fixed_column_counts)
    )
    # State is (added cells, number of odd final columns).
    dynamic = {(0, 0): 1}
    for fixed in columns:
        updated: dict[tuple[int, int], int] = {}
        for (added, odd), count in dynamic.items():
            remaining = degree - fixed_size - added
            for extra in range(min(order - fixed, remaining) + 1):
                key = (added + extra, odd + (fixed + extra) % 2)
                updated[key] = updated.get(key, 0) + count * comb(
                    order - fixed, extra
                )
        dynamic = updated
    return dynamic.get((degree - fixed_size, 1), 0)


def endpoint_quintic_incidence(selected: int) -> int:
    """Exact maximum incidence degree for degree-five, record-one supports."""

    return one_record_incidence(5, selected)


@lru_cache(maxsize=None)
def one_record_incidence_for_order(
    order: int, degree: int, selected: int
) -> int:
    """Maximum incidence for one-record supports at a specified order."""

    return max(
        one_record_extension_count(order, partition, degree)
        for partition in integer_partitions(selected)
        if len(partition) <= order and max(partition, default=0) <= order
    )


def one_record_incidence(degree: int, selected: int) -> int:
    """Maximum incidence for degree-``degree`` one-record supports."""

    return one_record_incidence_for_order(ORDER, degree, selected)


def endpoint_quintic_singleton_slice_energies(
    order: int,
) -> tuple[float, ...]:
    """Exact endpoint-quintic squared slices next to a singleton.

    The endpoint support has degree five and one odd label on the link to
    the singleton.  Entry ``k`` is the largest squared moment sum over
    endpoint supports containing a fixed ``k``-cell subset, with the
    singleton coordinate also optimized.
    """

    q = order
    if q < 4 or q & (q - 1):
        raise ValueError(
            ("degree-five slice requires power-of-two order at least four", q)
        )
    w0 = 1 / q**2
    w1 = 1 / (q**2 * (q - 1) ** 2)
    w2 = 4 / (q**2 * (q - 1) ** 2 * (q - 2) ** 2)

    zero_xor_four_sets = q * (q - 1) * (q - 2) // 24
    count_5 = q * comb(q, 5)
    count_41 = q**2 * (q - 1) * comb(q, 4)
    high_41 = q**2 * (q - 1) * zero_xor_four_sets
    count_32 = q * (q - 1) * comb(q, 3) * comb(q, 2)
    count_221 = q**2 * comb(q - 1, 2) * comb(q, 2) ** 2
    high_221 = (
        q**2 * comb(q - 1, 2) * (q - 1) * q**2 // 4
    )
    slice_0 = (
        (count_5 + high_41) * w0
        + (count_41 - high_41 + count_32 + high_221) * w1
        + (count_221 - high_221) * w2
    )
    slice_1 = 5 * slice_0 / q**2

    half_nonzero_xors = q // 2 - 1
    slice_2 = (
        comb(q - 2, 3) * w0
        + q
        * (q - 1)
        * (
            half_nonzero_xors * w0
            + (comb(q - 2, 2) - half_nonzero_xors) * w1
        )
        + (
            (q - 2) * (q - 1) * comb(q, 2)
            + (q - 1) * comb(q, 3)
        )
        * w1
        + (q - 1)
        * (q - 2)
        * q
        * (
            (q // 2) * w1
            + (comb(q, 2) - q // 2) * w2
        )
    )
    slice_3 = (
        comb(q - 3, 2) * w0
        + q * (q - 1) * (w0 + (q - 4) * w1)
        + (q - 1) * comb(q, 2) * w1
    )
    return (
        slice_0,
        slice_1,
        slice_2,
        slice_3,
        1 - 4 / q**2,
        1 / q**2,
    )


def middle_quintic_incidence_bound(
    order: int, selected: int
) -> int:
    """Incidence bound for degree-five supports with both records one.

    The zero- and one-cell entries use the safe one-record relaxation.
    The remaining entries are exact maxima.  Fixed pairs maximize in one
    row or column, fixed triples maximize on an L-shape, and fixed four
    cells maximize on a four-cycle.
    """

    q = order
    if q < 4:
        raise ValueError(("middle-quintic incidence requires order four", q))
    if selected in (0, 1):
        return one_record_incidence_for_order(q, 5, selected)
    if selected == 2:
        return (q - 2) * (q - 1) * (5 * q - 4)
    if selected == 3:
        return 2 * (q - 2) * (2 * q - 1)
    if selected == 4:
        return q**2 - 4
    if selected == 5:
        return 1
    raise ValueError(("invalid quintic split", selected))


def incidence_coefficient(
    profile: tuple[int, ...],
    split: tuple[int, ...],
    families: tuple,
    maximum_entry: float,
) -> float:
    row = prod(family(selected) for family, selected in zip(families, split, strict=True))
    column = prod(
        family(degree - selected)
        for family, degree, selected in zip(
            families, profile, split, strict=True
        )
    )
    return maximum_entry * min(sqrt(row), sqrt(column))


def record_three_endpoint_incidence(selected: int) -> int:
    """Triples with three distinct hidden labels and arbitrary partners."""

    q = ORDER
    values = (
        comb(q, 3) * q**3,
        comb(q - 1, 2) * q**2,
        q * (q - 2),
        1,
    )
    return values[selected]


def record_three_star_incidence(selected: int) -> int:
    """Triples with record sizes one on one side and three on the other."""

    q = ORDER
    values = (
        q * comb(q, 3) * (3 * q - 2),
        (q - 1) * (q - 2) * (3 * q - 2) // 2,
        q * (q - 2),
        1,
    )
    return values[selected]


def matching_cubic_incidence(selected: int) -> int:
    """Cubics with three distinct labels on both sides."""

    q = ORDER
    values = (
        6 * comb(q, 3) ** 2,
        2 * comb(q - 1, 2) ** 2,
        (q - 2) ** 2,
        1,
    )
    return values[selected]


def endpoint_singleton_slice_energies(order: int) -> tuple[float, ...]:
    """Exact squared slices for a cubic endpoint next to a singleton."""

    q = order
    return (
        (q * q + 2) / 6,
        (q * q + 2) / (2 * q * q),
        (q * q - 2 * q + 2) / (q * q * (q - 1)),
        1 / (q * q),
    )


def star_singleton_slice_energies(order: int) -> tuple[float, ...]:
    """Exact squared slices for a record-(1,3) cubic and singleton."""

    q = order
    return (
        (q * q - 4) / 6,
        (q * q - 4) / (2 * q * q),
        (q - 2) / (q * (q - 1)),
        1 / (q * q),
    )


def record_three_output_slice_energies(
    order: int,
) -> tuple[float, ...]:
    """Exact pure-record-three output slices for one fixed cubic.

    The output family contains every cubic with three distinct matching
    labels.  The four entries fix zero through three output cells.
    """

    q = order
    return (
        q * q / ((q - 1) * (q - 2)),
        3 / ((q - 1) * (q - 2)),
        12 / (q * (q - 1) ** 2 * (q - 2)),
        36 / (q * q * (q - 1) ** 2 * (q - 2) ** 2),
    )


def record_three_star_tail_slice_bounds(
    order: int,
) -> tuple[float, ...]:
    """Record-three slices after retaining the final singleton energy.

    A type-A star has singleton row energy one, a type-B star has row
    energy ``1/(q-1)^2``, and a matching cubic has zero singleton link.
    The pure record-three slice controls the type-B and matching remainder;
    maximum-entry incidence controls the exceptional type-A subfamily.
    """

    q = order
    pure = record_three_output_slice_energies(q)
    type_a_incidence = (
        q * comb(q, 3),
        comb(q - 1, 2),
        q - 2,
        1,
    )
    type_b_tail = 1 / (q - 1) ** 2
    maximum_record_three_entry = 1 / comb(q, 3)
    return tuple(
        type_b_tail * pure[selected]
        + (1 - type_b_tail)
        * maximum_record_three_entry**2
        * type_a_incidence[selected]
        for selected in range(4)
    )


def endpoint_to_star_output_slice_bounds(
    order: int,
) -> tuple[float, ...]:
    """Squared output-slice bounds for one fixed endpoint cubic."""

    q = order
    # Type A: all three matching labels coincide.  Type B: one label is
    # odd and a second label is repeated.  Every type-B entry is bounded by
    # the exceptional record-one value.
    type_a_incidence = (
        q * comb(q, 3),
        comb(q - 1, 2),
        q - 2,
        1,
    )
    type_b_incidence = (
        q * q * (q - 1) * comb(q - 1, 2),
        3 * (q - 1) ** 2 * (q - 2) / 2,
        (q - 1) * (q - 2),
        0,
    )
    exceptional_squared = (q + 2) ** 2 / (
        q * q * (q - 1) ** 2 * (q - 2) ** 2
    )
    return tuple(
        type_a_incidence[selected] / (q * q)
        + type_b_incidence[selected] * exceptional_squared
        if selected < 3
        else 1 / (q * q)
        for selected in range(4)
    )


def star_to_endpoint_input_slice_bounds(
    order: int,
) -> tuple[float, ...]:
    """Squared input-slice bounds for one fixed record-(1,3) cubic."""

    q = order
    endpoint_degrees = (
        q * comb(q, 3) + q * q * (q - 1) * comb(q, 2),
        comb(q - 1, 2)
        + (q - 1) * comb(q, 2)
        + q * (q - 1) ** 2,
        q * q - 2,
        1,
    )
    exceptional_squared = (q + 2) ** 2 / (
        q * q * (q - 1) ** 2 * (q - 2) ** 2
    )
    baseline_squared = 1 / (q * q * (q - 1) ** 2)
    exceptional_count = q**3 * (q - 2) ** 2 / 8
    type_b_full = (
        endpoint_degrees[0] * baseline_squared
        + exceptional_count * (exceptional_squared - baseline_squared)
    )
    type_b = tuple(
        min(type_b_full, exceptional_squared * degree)
        for degree in endpoint_degrees
    )
    type_a = endpoint_singleton_slice_energies(q)
    return tuple(max(a, b) for a, b in zip(type_a, type_b, strict=True))


def symmetric_row_column_coefficient(
    profile: tuple[int, ...],
    split: tuple[int, ...],
    row_energy,
) -> float:
    complement = tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )
    return sqrt(min(row_energy(split), row_energy(complement)))


def rank_incidence_coefficient(
    profile: tuple[int, ...],
    split: tuple[int, ...],
    families: tuple,
    maximum_entry: float,
) -> float:
    smaller_side = min(sum(split), sum(profile) - sum(split))
    return min(
        ORDER**smaller_side * maximum_entry,
        incidence_coefficient(profile, split, families, maximum_entry),
    )


def separated_triple_cubic_coefficient(
    split: tuple[int, ...],
) -> float:
    """The (3,1,3,3) record-one plus record-three sectors."""

    q = ORDER
    profile = (3, 1, 3, 3)
    record_one_entry = (q + 2) / (
        q**3 * (q - 1) ** 2 * (q - 2)
    )
    record_one = rank_incidence_coefficient(
        profile,
        split,
        (
            endpoint_record_one_incidence,
            singleton_incidence,
            l_shape_incidence,
            endpoint_record_one_incidence,
        ),
        record_one_entry,
    )
    endpoint = endpoint_singleton_slice_energies(q)
    star = star_singleton_slice_energies(q)
    record_three = record_three_output_slice_energies(q)

    def row_energy(local_split: tuple[int, ...]) -> float:
        # The singleton in block two is either fixed or summed over q^2
        # coordinates.  Conditional on it, the other three squared slices
        # separate exactly.
        singleton_factor = q ** (2 * (1 - local_split[1]))
        return (
            singleton_factor
            * endpoint[local_split[0]]
            * star[local_split[2]]
            * record_three[local_split[3]]
        )

    return record_one + symmetric_row_column_coefficient(
        profile, split, row_energy
    )


def end_chain_triple_cubic_coefficient(
    split: tuple[int, ...],
) -> float:
    """The four record triples of (3,3,3,1)."""

    q = ORDER
    profile = (3, 3, 3, 1)
    incidence_sectors = (
        (
            (
                endpoint_record_one_incidence,
                l_shape_incidence,
                l_shape_incidence,
                singleton_incidence,
            ),
            (q + 2) / (q**3 * (q - 1) ** 2 * (q - 2)),
        ),
        (
            (
                record_three_endpoint_incidence,
                record_three_star_incidence,
                l_shape_incidence,
                singleton_incidence,
            ),
            1 / (comb(q, 3) * q**2 * (q - 1)),
        ),
        (
            (
                record_three_endpoint_incidence,
                matching_cubic_incidence,
                record_three_star_incidence,
                singleton_incidence,
            ),
            1 / (q * comb(q, 3) ** 2),
        ),
    )
    coefficient = sum(
        rank_incidence_coefficient(
            profile, split, families, maximum_entry
        )
        for families, maximum_entry in incidence_sectors
    )

    endpoint_degrees = tuple(
        endpoint_record_one_incidence(selected) for selected in range(4)
    )
    star_degrees = tuple(
        record_three_star_incidence(selected) for selected in range(4)
    )
    endpoint_to_star = endpoint_to_star_output_slice_bounds(q)
    star_to_endpoint = star_to_endpoint_input_slice_bounds(q)
    joint_first_link = tuple(
        tuple(
            min(
                endpoint_degrees[left_selected]
                * endpoint_to_star[right_selected],
                star_degrees[right_selected]
                * star_to_endpoint[left_selected],
            )
            for right_selected in range(4)
        )
        for left_selected in range(4)
    )
    weighted_tail = record_three_star_tail_slice_bounds(q)

    def row_energy(local_split: tuple[int, ...]) -> float:
        return (
            joint_first_link[local_split[0]][local_split[1]]
            * weighted_tail[local_split[2]]
            / q ** (2 * local_split[3])
        )

    # This is the (1,3,1) record triple.  The final singleton energy is
    # already retained in ``weighted_tail`` rather than charged at its
    # worst type-A value for every support.
    coefficient += symmetric_row_column_coefficient(
        profile, split, row_energy
    )
    return coefficient


def triple_cubic_coefficient(
    profile: tuple[int, ...], split: tuple[int, ...]
) -> float:
    if profile == (3, 3, 3, 1):
        return end_chain_triple_cubic_coefficient(split)
    if profile == (1, 3, 3, 3):
        return end_chain_triple_cubic_coefficient(tuple(reversed(split)))
    if profile == (3, 1, 3, 3):
        return separated_triple_cubic_coefficient(split)
    if profile == (3, 3, 1, 3):
        return separated_triple_cubic_coefficient(tuple(reversed(split)))
    raise ValueError(("not a triple-cubic profile", profile))


def separated_quintic_cubic_coefficient(
    profile: tuple[int, ...], split: tuple[int, ...]
) -> float:
    """Fixed-split contraction for (5,1,3,1) and its reversal.

    The singleton neighbors force the cubic to be an L-shape and the
    quintic endpoint to have record one.  The three link entries are at
    most ``1/q``, ``1/[q(q-1)]``, and ``1/[q(q-1)]``.  Rank and exact
    support incidences are then combined cut by cut.
    """

    q = ORDER
    canonical = (5, 1, 3, 1)
    if profile == canonical:
        canonical_split = split
    elif profile == tuple(reversed(canonical)):
        canonical_split = tuple(reversed(split))
    else:
        raise ValueError(("not a separated quintic-cubic profile", profile))
    return rank_incidence_coefficient(
        canonical,
        canonical_split,
        (
            endpoint_quintic_incidence,
            singleton_incidence,
            l_shape_incidence,
            singleton_incidence,
        ),
        1 / (q**3 * (q - 1) ** 2),
    )


def separated_cubic_quintic_coefficient(
    profile: tuple[int, ...], split: tuple[int, ...]
) -> float:
    """Fixed-split contraction for (3,1,5,1) and its reversal.

    This bound uses the exact middle-quintic incidences for cuts fixing
    at least two quintic cells.  Its four-cycle cut is sharp for this local
    method, but the resulting coefficient still fits in the partial ledger.
    """

    q = ORDER
    canonical = (3, 1, 5, 1)
    if profile == canonical:
        canonical_split = split
    elif profile == tuple(reversed(canonical)):
        canonical_split = tuple(reversed(split))
    else:
        raise ValueError(("not a separated cubic-quintic profile", profile))
    return rank_incidence_coefficient(
        canonical,
        canonical_split,
        (
            endpoint_record_one_incidence,
            singleton_incidence,
            lambda selected: middle_quintic_incidence_bound(q, selected),
            singleton_incidence,
        ),
        1 / (q**3 * (q - 1) ** 2),
    )


def central_record_three_coefficient(split: tuple[int, ...]) -> float:
    """Entry/incidence row--column bound for (1,3,3,1), record three."""

    q = ORDER
    entry = 1 / (q * q * comb(q, 3))
    row_degrees = (
        singleton_incidence(split[0]),
        record_three_star_incidence(split[1]),
        record_three_star_incidence(split[2]),
        singleton_incidence(split[3]),
    )
    column_degrees = (
        singleton_incidence(1 - split[0]),
        record_three_star_incidence(3 - split[1]),
        record_three_star_incidence(3 - split[2]),
        singleton_incidence(1 - split[3]),
    )
    return entry * min(sqrt(prod(row_degrees)), sqrt(prod(column_degrees)))


def adjacent_record_three_coefficient(
    split: tuple[int, ...], endpoint_block: int, middle_block: int
) -> float:
    """Chained slice bound for an endpoint/middle record-three path.

    The result is the better of the inherited entry/incidence estimate and
    the exact singleton--star/output squared-slice factorization.
    """

    q = ORDER
    entry = 1 / (q * q * comb(q, 3))

    def degree(block: int, selected: int) -> int:
        if block == endpoint_block:
            return record_three_endpoint_incidence(selected)
        if block == middle_block:
            return record_three_star_incidence(selected)
        return singleton_incidence(selected)

    row_degrees = tuple(degree(block, split[block]) for block in range(4))
    column_degrees = tuple(
        degree(block, (3 if block in (endpoint_block, middle_block) else 1) - split[block])
        for block in range(4)
    )
    incidence = entry * min(
        sqrt(prod(row_degrees)),
        sqrt(prod(column_degrees)),
    )

    # Chained squared-slice refinement.  In canonical order the profile is
    # (singleton, singleton, record-(1,3) star, record-three endpoint).
    # For fixed star and endpoint partial supports, first sum the pure
    # record-three output using C_l, then sum the preceding singleton--star
    # link using T_k.  The M_11 square contributes N^(1-s_0-s_1).
    if (endpoint_block, middle_block) == (3, 2):
        canonical_split = split
    elif (endpoint_block, middle_block) == (0, 1):
        canonical_split = tuple(reversed(split))
    else:
        raise ValueError(("nonadjacent endpoint/middle blocks", split))
    singleton_selected = canonical_split[0] + canonical_split[1]
    star_selected = canonical_split[2]
    endpoint_selected = canonical_split[3]
    star = star_singleton_slice_energies(q)
    endpoint = record_three_output_slice_energies(q)
    row_energy = (
        q ** (2 * (1 - singleton_selected))
        * star[star_selected]
        * endpoint[endpoint_selected]
    )
    complement_energy = (
        q ** (2 * (singleton_selected - 1))
        * star[3 - star_selected]
        * endpoint[3 - endpoint_selected]
    )
    chained = sqrt(min(row_energy, complement_energy))
    return min(incidence, chained)


def degree_eight_coefficient(
    profile: tuple[int, ...],
    split: tuple[int, ...],
    include_record_three: bool = True,
) -> float:
    """Safe Round 2 fixed-split constants for every degree-eight profile."""

    q = ORDER
    selected_marks = sum(split)
    smaller_side = min(selected_marks, 8 - selected_marks)
    rank_factor = q**smaller_side
    decorated = tuple(
        block for block, degree in enumerate(profile) if degree > 1
    )
    if len(decorated) == 1:
        block = decorated[0]
        families = tuple(
            endpoint_quintic_incidence
            if current == block
            else singleton_incidence
            for current in range(4)
        )
        if block in (0, 3):
            # The full entry is at most q^-3.  Combine split-dependent rank
            # with the sharp Gram-dressed endpoint cap q^-1.
            return min(
                1 / q,
                rank_factor / q**3,
                incidence_coefficient(
                    profile, split, families, 1 / q**3
                ),
            )
        # Middle-quintic entries are at most 1/[q^3(q-1)^2].
        entry = 1 / (q**3 * (q - 1) ** 2)
        return min(
            rank_factor * entry,
            incidence_coefficient(profile, split, families, entry),
        )

    if profile == DOUBLE_ENDPOINT_PROFILE:
        return double_endpoint_coefficient(split)
    if decorated == (0, 1):
        record_one_entry = (q + 2) / (
            q**3 * (q - 1) ** 2 * (q - 2)
        )
        record_one = min(
            adjacent_record_one_coefficient(split[0], split[1]),
            rank_factor * record_one_entry,
            incidence_coefficient(
                profile,
                split,
                (
                    endpoint_record_one_incidence,
                    l_shape_incidence,
                    singleton_incidence,
                    singleton_incidence,
                ),
                record_one_entry,
            ),
        )
        return (
            record_one
            + (
                adjacent_record_three_coefficient(split, 0, 1)
                if include_record_three
                else 0.0
            )
        )
    if decorated == (2, 3):
        record_one_entry = (q + 2) / (
            q**3 * (q - 1) ** 2 * (q - 2)
        )
        record_one = min(
            adjacent_record_one_coefficient(split[3], split[2]),
            rank_factor * record_one_entry,
            incidence_coefficient(
                profile,
                split,
                (
                    singleton_incidence,
                    singleton_incidence,
                    l_shape_incidence,
                    endpoint_record_one_incidence,
                ),
                record_one_entry,
            ),
        )
        return (
            record_one
            + (
                adjacent_record_three_coefficient(split, 3, 2)
                if include_record_three
                else 0.0
            )
        )
    if decorated == (1, 2):
        # Triangle over the central record-one and record-three sectors.
        record_one_entry = 1 / (q**3 * (q - 1) ** 2)
        record_one = min(
            rank_factor * record_one_entry,
            incidence_coefficient(
                profile,
                split,
                (
                    singleton_incidence,
                    l_shape_incidence,
                    l_shape_incidence,
                    singleton_incidence,
                ),
                record_one_entry,
            ),
        )
        record_three = (
            central_record_three_coefficient(split)
            if include_record_three
            else 0.0
        )
        return record_one + record_three
    # The two separated cubic profiles have a forced middle L-shape.
    record_one_entry = 1 / (q**3 * (q - 1) ** 2)
    if decorated == (0, 2):
        families = (
            endpoint_record_one_incidence,
            singleton_incidence,
            l_shape_incidence,
            singleton_incidence,
        )
    else:
        families = (
            singleton_incidence,
            l_shape_incidence,
            singleton_incidence,
            endpoint_record_one_incidence,
        )
    return min(
        rank_factor * record_one_entry,
        incidence_coefficient(
            profile, split, families, record_one_entry
        ),
    )


def multiplicity(state: tuple[int, ...], selected: tuple[int, ...]) -> int:
    return prod(
        comb(occupation, count)
        for occupation, count in zip(state, selected, strict=True)
    )


def paired_state(
    state: tuple[int, ...],
    profile: tuple[int, ...],
    split: tuple[int, ...],
) -> tuple[int, ...]:
    return tuple(
        occupation + degree - 2 * selected
        for occupation, degree, selected in zip(
            state, profile, split, strict=True
        )
    )


def certificate(
    include_degree_eight: bool = True,
    include_record_three: bool = True,
    include_known_high_degree: bool = True,
    beta: float = BETA,
    high_degree_coefficient: float | None = None,
    profile_split_coefficients: Mapping[
        tuple[tuple[int, ...], tuple[int, ...]], float
    ]
    | None = None,
) -> CompatibleOccupationResult:
    states = occupation_states()
    state_index = {state: index for index, state in enumerate(states)}
    profiles = (MINIMAL_PROFILE,) + DEGREE_SIX_PROFILES
    if include_degree_eight:
        profiles += DEGREE_EIGHT_PROFILES
    if include_known_high_degree:
        profiles += KNOWN_HIGH_DEGREE_PROFILES
    if high_degree_coefficient is not None:
        profiles += tuple(
            profile
            for profile in HIGH_DEGREE_PROFILES
            if profile not in profiles
        )
    mapped_profiles = (
        tuple(
            sorted(
                {profile for profile, _ in profile_split_coefficients}
            )
        )
        if profile_split_coefficients is not None
        else ()
    )
    profiles += tuple(
        profile for profile in mapped_profiles if profile not in profiles
    )
    term_left: list[int] = []
    term_right: list[int] = []
    term_constants: list[float] = []
    term_profiles: list[tuple[int, ...]] = []
    term_splits: list[tuple[int, ...]] = []
    term_unattenuated: list[float] = []

    for profile in profiles:
        attenuation = beta ** sum(profile)
        for split in profile_splits(profile):
            complement = tuple(
                degree - selected
                for degree, selected in zip(profile, split, strict=True)
            )
            if profile in mapped_profiles:
                unattenuated = float(
                    profile_split_coefficients.get((profile, split), 0.0)
                )
                if unattenuated < 0:
                    raise ValueError(
                        ("negative profile-split coefficient", profile, split)
                    )
            else:
                unattenuated = coefficient(
                    profile,
                    split,
                    include_record_three,
                    high_degree_coefficient,
                )
            local = unattenuated * attenuation
            for left_index, state in enumerate(states):
                if any(
                    occupation < selected
                    for occupation, selected in zip(state, split, strict=True)
                ):
                    continue
                partner = paired_state(state, profile, split)
                right_index = state_index.get(partner)
                if right_index is None:
                    continue
                left_count = multiplicity(state, split)
                right_count = multiplicity(partner, complement)
                if not left_count or not right_count:
                    continue
                term_left.append(left_index)
                term_right.append(right_index)
                term_constants.append(local * sqrt(left_count * right_count))
                term_profiles.append(profile)
                term_splits.append(split)
                term_unattenuated.append(unattenuated)

    # Aggregate all terms that join the same two occupation states.  With
    # x_i=sqrt(rho_i), the objective is exactly the Rayleigh quotient of a
    # nonnegative symmetric 210-by-210 matrix.  Its global maximum is the
    # Perron eigenvalue; no nonlinear solver is needed.
    edge_constants: dict[tuple[int, int], float] = {}
    for left, right, constant in zip(
        term_left, term_right, term_constants, strict=True
    ):
        edge = (min(left, right), max(left, right))
        edge_constants[edge] = edge_constants.get(edge, 0.0) + constant
    matrix = np.zeros((len(states), len(states)))
    for (left, right), constant in edge_constants.items():
        if left == right:
            raise AssertionError(("odd profile produced a loop", states[left]))
        matrix[left, right] += constant / 2
        matrix[right, left] += constant / 2

    # A tiny positive regularization selects a strictly positive Perron
    # vector even if the occupation graph is reducible.  Collatz--Wielandt
    # below is applied to the unregularized matrix and remains a valid upper.
    regularization = 1e-14
    _, eigenvectors = np.linalg.eigh(
        matrix + regularization * np.ones_like(matrix)
    )
    perron = np.abs(eigenvectors[:, -1])
    perron /= np.linalg.norm(perron)
    candidate = perron * perron
    left_values = candidate[np.asarray(term_left)]
    right_values = candidate[np.asarray(term_right)]
    constants = np.asarray(term_constants)
    values = constants * np.sqrt(left_values * right_values)
    direct_objective = float(values.sum())
    rayleigh = float(perron @ matrix @ perron)
    if abs(rayleigh - direct_objective) > 3e-10:
        raise AssertionError((rayleigh, direct_objective))

    gradient = (matrix @ perron) / perron
    supporting_upper = float(np.max(gradient) + 1e-11)
    if direct_objective > supporting_upper:
        raise AssertionError((direct_objective, supporting_upper))

    profile_totals: dict[tuple[int, ...], float] = {}
    for profile, value in zip(term_profiles, values, strict=True):
        profile_totals[profile] = profile_totals.get(profile, 0.0) + float(
            value
        )
    return CompatibleOccupationResult(
        objective=direct_objective,
        supporting_upper=supporting_upper,
        profile_contributions=tuple(
            sorted(
                profile_totals.items(),
                key=lambda item: item[1],
                reverse=True,
            )
        ),
        occupation_weights=tuple(
            (state, float(weight))
            for state, weight in zip(states, candidate, strict=True)
        ),
        support=tuple(
            sorted(
                (
                    (float(weight), state)
                    for state, weight in zip(states, candidate, strict=True)
                    if weight > 2e-7
                ),
                reverse=True,
            )
        ),
        maximum_gradient_state=states[int(np.argmax(gradient))],
        leading_terms=tuple(
            sorted(
                (
                    (
                        float(value),
                        profile,
                        split,
                        states[left],
                        states[right],
                        unattenuated,
                    )
                    for value, profile, split, left, right, unattenuated in zip(
                        values,
                        term_profiles,
                        term_splits,
                        term_left,
                        term_right,
                        term_unattenuated,
                        strict=True,
                    )
                ),
                reverse=True,
            )[:20]
        ),
    )


def main() -> None:
    result = certificate()
    promise = promise_concentration(ORDER * ORDER, BETA)
    partial_total = result.supporting_upper + promise.two_hypothesis_loss
    print(
        "occupation-compatible sector certificate: "
        f"objective={result.objective:.15g},"
        f"supporting_upper={result.supporting_upper:.15g},"
        f"beta={BETA:.15g},"
        f"promise_loss={promise.two_hypothesis_loss:.15g},"
        f"partial_total={partial_total:.15g},"
        f"threshold_slack={TV_THRESHOLD-partial_total:.15g},"
        f"maximum_gradient_state={result.maximum_gradient_state}"
    )
    print(
        "profile_contributions="
        + ";".join(
            f"{profile}:{value:.12g}"
            for profile, value in result.profile_contributions
        )
    )
    print(
        "support="
        + ";".join(
            f"{state}:{weight:.12g}" for weight, state in result.support
        )
    )
    print(
        "leading_terms="
        + ";".join(
            f"{profile}/{split}/{left}->{right}:"
            f"{value:.12g}@{unattenuated:.12g}"
            for value, profile, split, left, right, unattenuated in (
                result.leading_terms
            )
        )
    )


if __name__ == "__main__":
    main()
