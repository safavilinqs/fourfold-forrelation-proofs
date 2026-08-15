#!/usr/bin/env python3
"""Regression for compatible occupations and triple-cubic chained slices."""

from __future__ import annotations

from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from attenuation_promise_concentration import promise_concentration
from occupation_compatible_sector_optimization import (
    BETA,
    SEPARATED_CUBIC_QUINTIC_PROFILES,
    SEPARATED_QUINTIC_CUBIC_PROFILES,
    TRIPLE_CUBIC_PROFILES,
    adjacent_record_three_coefficient,
    certificate,
    endpoint_singleton_slice_energies,
    endpoint_quintic_singleton_slice_energies,
    endpoint_to_star_output_slice_bounds,
    integer_partitions,
    middle_quintic_incidence_bound,
    one_record_extension_count,
    profile_splits,
    record_three_output_slice_energies,
    record_three_star_tail_slice_bounds,
    separated_cubic_quintic_coefficient,
    separated_quintic_cubic_coefficient,
    star_singleton_slice_energies,
    star_to_endpoint_input_slice_bounds,
    triple_cubic_coefficient,
)
from opposite_endpoint_orbit_scan import (
    aligned_vertical_orbit_coefficient,
    endpoint_moment,
    exhaustive_fixed_orbit_maximum,
    orbit_block,
)


def record_size(support: tuple[int, ...], order: int, axis: int) -> int:
    counts = Counter(divmod(coordinate, order)[axis] for coordinate in support)
    return sum(count % 2 for count in counts.values())


def incidence_degrees(
    supports: list[tuple[int, ...]], dimension: int, degree: int
) -> tuple[int, ...]:
    result = []
    for selected in range(degree + 1):
        counts: Counter[tuple[int, ...]] = Counter()
        for support in supports:
            counts.update(combinations(support, selected))
        result.append(max(counts.values()))
    return tuple(result)


def incidence_matrix(
    supports: list[tuple[int, ...]], dimension: int, selected: int
) -> np.ndarray:
    parts = list(combinations(range(dimension), selected))
    part_index = {part: index for index, part in enumerate(parts)}
    result = np.zeros((len(parts), len(supports)))
    for column, support in enumerate(supports):
        for part in combinations(support, selected):
            result[part_index[part], column] = 1
    return result


def both_record_extension_count(
    order: int, fixed: tuple[int, ...], degree: int = 5
) -> int:
    """Exact small-order extension count for fixed bipartite edges."""

    available = tuple(
        coordinate
        for coordinate in range(order * order)
        if coordinate not in fixed
    )
    return sum(
        record_size(fixed + extra, order, 0) == 1
        and record_size(fixed + extra, order, 1) == 1
        for extra in combinations(available, degree - len(fixed))
    )


def q4_link_moments() -> tuple[
    list[tuple[int, ...]],
    list[tuple[int, ...]],
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    order = 4
    hadamard = np.array([[1]], dtype=np.int8)
    while hadamard.shape[0] < order:
        hadamard = np.block(
            [[hadamard, hadamard], [hadamard, -hadamard]]
        )
    left_values = []
    right_values = []
    for permutation in permutations(range(order)):
        for signs in product((-1, 1), repeat=order):
            signed_permutation = np.zeros((order, order), dtype=np.int8)
            for column, row in enumerate(permutation):
                signed_permutation[row, column] = signs[column]
            left_values.append((hadamard @ signed_permutation).reshape(-1))
            right_values.append((signed_permutation @ hadamard).reshape(-1))
    left = np.asarray(left_values, dtype=np.int8)
    right = np.asarray(right_values, dtype=np.int8)
    degree_one = list(combinations(range(order * order), 1))
    degree_three = list(combinations(range(order * order), 3))
    degree_five = list(combinations(range(order * order), 5))

    def features(
        values: np.ndarray, supports: list[tuple[int, ...]]
    ) -> np.ndarray:
        return np.asarray(
            [np.prod(values[:, support], axis=1) for support in supports]
        ).T

    left_one = features(left, degree_one)
    right_one = features(right, degree_one)
    left_three = features(left, degree_three)
    right_three = features(right, degree_three)
    left_five = features(left, degree_five)
    right_five = features(right, degree_five)
    normalization = len(left)
    m13 = left_one.T.astype(float) @ right_three / normalization
    m31 = left_three.T.astype(float) @ right_one / normalization
    m33 = left_three.T.astype(float) @ right_three / normalization
    m15 = left_one.T.astype(float) @ right_five / normalization
    m51 = left_five.T.astype(float) @ right_one / normalization
    return degree_three, degree_five, m13, m31, m33, m15, m51


def maximum_weighted_slice(
    supports: list[tuple[int, ...]],
    weights: np.ndarray,
    selected: int,
) -> float:
    """Largest weighted extension sum without a dense incidence matrix."""

    totals: dict[tuple[int, ...], np.ndarray] = {}
    for support, row in zip(supports, weights, strict=True):
        for part in combinations(support, selected):
            if part not in totals:
                totals[part] = np.array(row, copy=True)
            else:
                totals[part] += row
    return max(float(np.max(total)) for total in totals.values())


def q4_incidence_checks() -> None:
    order = 4
    dimension = order * order
    for degree in (3, 5):
        supports = [
            support
            for support in combinations(range(dimension), degree)
            if record_size(support, order, 1) == 1
        ]
        observed = incidence_degrees(supports, dimension, degree)
        expected = tuple(
            max(
                one_record_extension_count(order, partition, degree)
                for partition in integer_partitions(selected)
                if len(partition) <= order
            )
            for selected in range(degree + 1)
        )
        if observed != expected:
            raise AssertionError(("one-record incidence", degree, observed, expected))

    cubics = list(combinations(range(dimension), 3))
    endpoint_one = [
        support
        for support in cubics
        if record_size(support, order, 1) == 1
    ]
    l_shapes = [
        support
        for support in endpoint_one
        if record_size(support, order, 0) == 1
    ]
    endpoint_three = [
        support
        for support in cubics
        if record_size(support, order, 1) == 3
    ]
    stars = [
        support
        for support in endpoint_three
        if record_size(support, order, 0) == 1
    ]
    matchings = [
        support
        for support in endpoint_three
        if record_size(support, order, 0) == 3
    ]
    expected = {
        "endpoint_one": (304, 57, 14, 1),
        "l_shape": (144, 27, 6, 1),
        "endpoint_three": (256, 48, 8, 1),
        "star": (160, 30, 8, 1),
        "matching": (96, 18, 4, 1),
    }
    observed = {
        "endpoint_one": incidence_degrees(endpoint_one, dimension, 3),
        "l_shape": incidence_degrees(l_shapes, dimension, 3),
        "endpoint_three": incidence_degrees(endpoint_three, dimension, 3),
        "star": incidence_degrees(stars, dimension, 3),
        "matching": incidence_degrees(matchings, dimension, 3),
    }
    if observed != expected:
        raise AssertionError(("cubic incidence", observed, expected))


def middle_quintic_formula_checks() -> None:
    """Check every fixed-pair/triple isomorphism formula exactly."""

    for order in range(4, 9):
        pair_counts = {
            "same_row": both_record_extension_count(order, (0, 1)),
            "matching": both_record_extension_count(
                order, (0, order + 1)
            ),
        }
        expected_pairs = {
            "same_row": (order - 2) * (order - 1) * (5 * order - 4),
            "matching": (order - 2) * (15 * order - 22),
        }
        if pair_counts != expected_pairs:
            raise AssertionError(
                ("middle-quintic pair formulas", order, pair_counts)
            )

        triple_counts = {
            "star": both_record_extension_count(order, (0, 1, 2)),
            "L": both_record_extension_count(
                order, (0, 1, order)
            ),
            "star_disjoint": both_record_extension_count(
                order, (0, 1, 2 * order + 2)
            ),
            "matching": both_record_extension_count(
                order, (0, order + 1, 2 * order + 2)
            ),
        }
        expected_triples = {
            "star": 3 * (order - 1),
            "L": 2 * (order - 2) * (2 * order - 1),
            "star_disjoint": 7 * order - 11,
            "matching": 9,
        }
        if triple_counts != expected_triples:
            raise AssertionError(
                ("middle-quintic triple formulas", order, triple_counts)
            )


def q4_slice_checks() -> None:
    order = 4
    dimension = order * order
    (
        degree_three,
        degree_five,
        m13,
        m31,
        m33,
        m15,
        m51,
    ) = q4_link_moments()
    formula_m31 = np.asarray(
        [
            [endpoint_moment(support, z, order, 3, False) for z in range(dimension)]
            for support in degree_three
        ]
    )
    formula_m51 = np.asarray(
        [
            [endpoint_moment(support, z, order, 5, False) for z in range(dimension)]
            for support in degree_five
        ]
    )
    if not np.array_equal(formula_m31, m31) or not np.array_equal(
        formula_m51, m51
    ):
        raise AssertionError("xor-labelled endpoint formulas")
    row_record = np.asarray(
        [record_size(support, order, 0) for support in degree_three]
    )
    column_record = np.asarray(
        [record_size(support, order, 1) for support in degree_three]
    )
    endpoint_mask = column_record == 1
    star_mask = (row_record == 1) & (column_record == 3)
    record_three_mask = row_record == 3
    endpoint_supports = [
        support
        for support, keep in zip(
            degree_three, endpoint_mask, strict=True
        )
        if keep
    ]
    star_supports = [
        support
        for support, keep in zip(degree_three, star_mask, strict=True)
        if keep
    ]
    record_three_supports = [
        support
        for support, keep in zip(
            degree_three, record_three_mask, strict=True
        )
        if keep
    ]
    endpoint_incidence = tuple(
        incidence_matrix(endpoint_supports, dimension, selected)
        for selected in range(4)
    )
    star_incidence = tuple(
        incidence_matrix(star_supports, dimension, selected)
        for selected in range(4)
    )
    record_three_incidence = tuple(
        incidence_matrix(record_three_supports, dimension, selected)
        for selected in range(4)
    )

    star_singleton_squared = np.abs(m13[:, star_mask]) ** 2
    observed_star = tuple(
        float(np.max(star_singleton_squared @ incidence.T))
        for incidence in star_incidence
    )
    expected_star = star_singleton_slice_energies(order)
    if not np.allclose(observed_star, expected_star, atol=2e-12):
        raise AssertionError(("star singleton slices", observed_star, expected_star))

    record_three_squared = np.abs(
        m33[np.ix_(star_mask, record_three_mask)]
    ) ** 2
    observed_record_three = tuple(
        float(np.max(record_three_squared @ incidence.T))
        for incidence in record_three_incidence
    )
    expected_record_three = record_three_output_slice_energies(order)
    if not np.allclose(
        observed_record_three, expected_record_three, atol=2e-12
    ):
        raise AssertionError(
            (
                "record-three output slices",
                observed_record_three,
                expected_record_three,
            )
        )

    # Compose the two exact squared-slice tables on the dominant adjacent
    # double-cubic cut.  This directly checks the chained inequality at q=4
    # before it is used at q=32.
    for star_selected, endpoint_selected in ((2, 1), (1, 2)):
        endpoint_tail = (
            record_three_squared
            @ record_three_incidence[endpoint_selected].T
        )
        observed_joint = 0.0
        for singleton in range(dimension):
            weighted_star = (
                star_incidence[star_selected]
                * star_singleton_squared[singleton][None, :]
            )
            observed_joint = max(
                observed_joint,
                float(np.max(weighted_star @ endpoint_tail)),
            )
        chained_upper = (
            expected_star[star_selected]
            * expected_record_three[endpoint_selected]
        )
        if observed_joint > chained_upper + 2e-12:
            raise AssertionError(
                (
                    "adjacent record-three chained slice",
                    star_selected,
                    endpoint_selected,
                    observed_joint,
                    chained_upper,
                )
            )

    q32_adjacent = adjacent_record_three_coefficient(
        (0, 1, 2, 1),
        endpoint_block=3,
        middle_block=2,
    )
    if not np.isclose(q32_adjacent, 0.002545238139708185, atol=2e-14):
        raise AssertionError(("q32 adjacent chained coefficient", q32_adjacent))
    reversed_adjacent = adjacent_record_three_coefficient(
        (1, 2, 1, 0),
        endpoint_block=0,
        middle_block=1,
    )
    if not np.isclose(q32_adjacent, reversed_adjacent, atol=2e-14):
        raise AssertionError(
            ("adjacent chained reversal", q32_adjacent, reversed_adjacent)
        )

    # Retain the actual final singleton row energy of every record-three
    # output support.  Matching supports have zero, type-B stars have
    # 1/(q-1)^2, and type-A stars have one.
    singleton_tail = np.sum(np.abs(m31[record_three_mask, :]) ** 2, axis=1)
    weighted = record_three_squared * singleton_tail[None, :]
    observed_weighted = tuple(
        float(np.max(weighted @ incidence.T))
        for incidence in record_three_incidence
    )
    expected_weighted_upper = record_three_star_tail_slice_bounds(order)
    if any(
        observed > expected + 2e-12
        for observed, expected in zip(
            observed_weighted, expected_weighted_upper, strict=True
        )
    ):
        raise AssertionError(
            (
                "record-three star-tail slices",
                observed_weighted,
                expected_weighted_upper,
            )
        )

    endpoint_to_star_squared = np.abs(
        m33[np.ix_(endpoint_mask, star_mask)]
    ) ** 2
    observed_output = tuple(
        float(np.max(endpoint_to_star_squared @ incidence.T))
        for incidence in star_incidence
    )
    output_upper = endpoint_to_star_output_slice_bounds(order)
    if any(
        observed > expected + 2e-12
        for observed, expected in zip(
            observed_output, output_upper, strict=True
        )
    ):
        raise AssertionError(
            ("endpoint-to-star output slices", observed_output, output_upper)
        )

    observed_input = tuple(
        float(np.max(incidence @ endpoint_to_star_squared))
        for incidence in endpoint_incidence
    )
    input_upper = star_to_endpoint_input_slice_bounds(order)
    if any(
        observed > expected + 2e-12
        for observed, expected in zip(
            observed_input, input_upper, strict=True
        )
    ):
        raise AssertionError(
            ("star-to-endpoint input slices", observed_input, input_upper)
        )

    row_record_five = np.asarray(
        [record_size(support, order, 0) for support in degree_five]
    )
    column_record_five = np.asarray(
        [record_size(support, order, 1) for support in degree_five]
    )
    endpoint_five_mask = column_record_five == 1
    endpoint_five_supports = [
        support
        for support, keep in zip(
            degree_five, endpoint_five_mask, strict=True
        )
        if keep
    ]
    endpoint_five_squared = np.abs(m51[endpoint_five_mask, :]) ** 2
    observed_quintic = tuple(
        maximum_weighted_slice(
            endpoint_five_supports, endpoint_five_squared, selected
        )
        for selected in range(6)
    )
    expected_quintic = endpoint_quintic_singleton_slice_energies(order)
    if not np.allclose(observed_quintic, expected_quintic, atol=2e-12):
        raise AssertionError(
            ("endpoint-quintic singleton slices", observed_quintic, expected_quintic)
        )

    middle_five_mask = (row_record_five == 1) & (column_record_five == 1)
    middle_five_supports = [
        support
        for support, keep in zip(
            degree_five, middle_five_mask, strict=True
        )
        if keep
    ]
    observed_middle_incidence = incidence_degrees(
        middle_five_supports, dimension, 5
    )
    expected_middle_incidence = (1008, 315, 96, 28, 12, 1)
    if observed_middle_incidence != expected_middle_incidence:
        raise AssertionError(
            (
                "middle-quintic incidence",
                observed_middle_incidence,
                expected_middle_incidence,
            )
        )
    incidence_bounds = tuple(
        middle_quintic_incidence_bound(order, selected)
        for selected in range(6)
    )
    if any(
        observed > bound
        for observed, bound in zip(
            observed_middle_incidence, incidence_bounds, strict=True
        )
    ) or observed_middle_incidence[2:] != incidence_bounds[2:]:
        raise AssertionError(
            (
                "middle-quintic incidence bounds",
                observed_middle_incidence,
                incidence_bounds,
            )
        )

    # A fixed four-cycle has q^2-4 compatible fifth-edge completions.  Every
    # completion simultaneously saturates the two adjacent singleton-link
    # entry bounds, so this D4 incidence cannot be improved locally.
    cycle = frozenset((0, 1, order, order + 1))
    cycle_indices = [
        index
        for index, support in enumerate(degree_five)
        if cycle.issubset(support) and middle_five_mask[index]
    ]
    if len(cycle_indices) != order**2 - 4:
        raise AssertionError(("middle-quintic four-cycle", len(cycle_indices)))
    link_bound = 1 / (order * (order - 1))
    for index in cycle_indices:
        if not np.isclose(
            np.max(np.abs(m15[:, index])), link_bound, atol=2e-12
        ) or not np.isclose(
            np.max(np.abs(m51[index, :])), link_bound, atol=2e-12
        ):
            raise AssertionError(("middle-quintic saturated link", index))


def opposite_endpoint_orbit_checks() -> None:
    for order in (4, 8, 16, 32):
        for high_only in (True, False):
            observed = orbit_block(
                order,
                order,
                order,
                (0, order, 1),
                high_only=high_only,
            ).normalized_nuclear
            expected = aligned_vertical_orbit_coefficient(order, high_only)
            if not np.isclose(observed, expected, atol=2e-14):
                raise AssertionError(
                    ("opposite endpoint aligned orbit", order, observed, expected)
                )

    # At q=4 all translation orbits can be scanned exactly.  The closed
    # value is the largest fixed-orbit value in both endpoint calculations.
    order = 4
    for high_only in (True, False):
        observed = exhaustive_fixed_orbit_maximum(
            order, high_only
        ).coefficient
        expected = aligned_vertical_orbit_coefficient(order, high_only)
        if not np.isclose(observed, expected, atol=2e-14):
            raise AssertionError(
                ("q4 fixed-orbit exhaustive maximum", observed, expected)
            )


def triple_cubic_symmetry_checks() -> None:
    for profile in ((3, 3, 3, 1), (3, 1, 3, 3)):
        reverse = tuple(reversed(profile))
        for split in profile_splits(profile):
            left = triple_cubic_coefficient(profile, split)
            right = triple_cubic_coefficient(
                reverse, tuple(reversed(split))
            )
            if not np.isclose(left, right, atol=2e-14):
                raise AssertionError(
                    ("triple-cubic reversal", profile, split, left, right)
                )

    maxima = {
        profile: max(
            triple_cubic_coefficient(profile, split)
            for split in profile_splits(profile)
        )
        for profile in TRIPLE_CUBIC_PROFILES
    }
    if not 0.03428 < maxima[(3, 3, 3, 1)] < 0.03429:
        raise AssertionError(("end-chain maximum", maxima))
    if not 0.05419 < maxima[(3, 1, 3, 3)] < 0.05420:
        raise AssertionError(("separated maximum", maxima))


def quintic_cubic_checks() -> None:
    for profile in SEPARATED_QUINTIC_CUBIC_PROFILES:
        reverse = tuple(reversed(profile))
        for split in profile_splits(profile):
            left = separated_quintic_cubic_coefficient(profile, split)
            right = separated_quintic_cubic_coefficient(
                reverse, tuple(reversed(split))
            )
            if not np.isclose(left, right, atol=2e-14):
                raise AssertionError(
                    ("quintic-cubic reversal", profile, split, left, right)
                )
    for profile in SEPARATED_CUBIC_QUINTIC_PROFILES:
        reverse = tuple(reversed(profile))
        for split in profile_splits(profile):
            left = separated_cubic_quintic_coefficient(profile, split)
            right = separated_cubic_quintic_coefficient(
                reverse, tuple(reversed(split))
            )
            if not np.isclose(left, right, atol=2e-14):
                raise AssertionError(
                    ("cubic-quintic reversal", profile, split, left, right)
                )
    maximum = max(
        separated_quintic_cubic_coefficient((5, 1, 3, 1), split)
        for split in profile_splits((5, 1, 3, 1))
    )
    if not 0.02467 < maximum < 0.02468:
        raise AssertionError(("quintic-cubic maximum", maximum))

    middle_maximum = max(
        separated_cubic_quintic_coefficient((3, 1, 5, 1), split)
        for split in profile_splits((3, 1, 5, 1))
    )
    if not 0.03320 < middle_maximum < 0.03321:
        raise AssertionError(
            ("middle-quintic incidence maximum", middle_maximum)
        )

    q32_slices = endpoint_quintic_singleton_slice_energies(32)
    expected_q32 = (
        52685,
        263425 / 1024,
        159457 / 7936,
        22365 / 15872,
        255 / 256,
        1 / 1024,
    )
    if not np.allclose(q32_slices, expected_q32, atol=2e-12):
        raise AssertionError(("q32 endpoint-quintic slices", q32_slices))
    if tuple(
        middle_quintic_incidence_bound(32, selected)
        for selected in range(2, 6)
    ) != (145080, 3780, 1020, 1):
        raise AssertionError("q32 middle-quintic incidence")

    # Directly pairing the cubic and quintic endpoint slice energies in the
    # (3,1,1,5) chain is far too weak.  This permanent negative result forces
    # a compound contraction rather than another row-energy application.
    cubic_slices = endpoint_singleton_slice_energies(32)
    naive_double_endpoint = max(
        (
            min(
                cubic_slices[cubic_selected]
                * q32_slices[quintic_selected],
                cubic_slices[3 - cubic_selected]
                * q32_slices[5 - quintic_selected],
            )
        )
        ** 0.5
        for cubic_selected in range(4)
        for quintic_selected in range(6)
    )
    if not 0.78032 < naive_double_endpoint < 0.78034:
        raise AssertionError(
            ("cubic-quintic row-energy barrier", naive_double_endpoint)
        )


def occupation_pairing_check() -> None:
    dimension = 5
    supports = [
        frozenset(chosen)
        for size in range(4)
        for chosen in combinations(range(dimension), size)
    ]
    for left in supports:
        for right in supports:
            profile = len(left ^ right)
            selected = len(left - right)
            if len(left) - selected != len(right) - (profile - selected):
                raise AssertionError((left, right, profile, selected))


def main() -> None:
    q4_incidence_checks()
    middle_quintic_formula_checks()
    q4_slice_checks()
    opposite_endpoint_orbit_checks()
    triple_cubic_symmetry_checks()
    quintic_cubic_checks()
    occupation_pairing_check()

    degree_eight_beta = 313 / 400
    result = certificate(
        beta=degree_eight_beta, include_known_high_degree=False
    )
    if not 0.23248 < result.objective < result.supporting_upper < 0.23249:
        raise AssertionError(("degree-eight compatible ledger", result))
    promise = promise_concentration(1024, degree_eight_beta)
    total = result.supporting_upper + promise.two_hypothesis_loss
    if not total < 0.249 < 1 / 3:
        raise AssertionError(("partial finite-size budget", total, result, promise))

    # The proved high endpoints, four triple cubics, and both separated
    # cubic--quintic reversal pairs now fit with every sector through degree
    # eight.  Other degree-ten/twelve profiles remain omitted.
    known = certificate(beta=BETA)
    known_promise = promise_concentration(1024, BETA)
    known_total = known.supporting_upper + known_promise.two_hypothesis_loss
    if not 0.2797 < known_total < 0.2798 < 1 / 3:
        raise AssertionError(("known high-degree partial ledger", known_total, known))

    # After charging both separated cubic--quintic pairs at their proved
    # coefficients, a common coefficient 1/32 for every still-open
    # degree-10/12 profile would pass, while 1/24 would not.  These are
    # target diagnostics, not claims that either coefficient is proved.
    passing_target = certificate(
        beta=BETA, high_degree_coefficient=1 / 32
    )
    failing_target = certificate(
        beta=BETA, high_degree_coefficient=1 / 24
    )
    passing_total = (
        passing_target.supporting_upper
        + known_promise.two_hypothesis_loss
    )
    failing_total = (
        failing_target.supporting_upper
        + known_promise.two_hypothesis_loss
    )
    if not passing_total < 1 / 3 < failing_total:
        raise AssertionError(("high-degree target", passing_total, failing_total))

    print(
        "occupation-compatible sector ledger passed: "
        f"degree8_upper={result.supporting_upper:.12g},"
        f"promise_loss={promise.two_hypothesis_loss:.12g},"
        f"partial_total={total:.12g},"
        f"partial_slack={1/3-total:.12g},"
        f"known_high_total={known_total:.12g},"
        f"known_high_slack={1/3-known_total:.12g},"
        f"target_1_over_32={passing_total:.12g},"
        f"target_1_over_24={failing_total:.12g}"
    )


if __name__ == "__main__":
    main()
