#!/usr/bin/env python3
"""Regression for the q64 noncubic and recovered-universal insertion."""

from __future__ import annotations

from itertools import combinations, permutations, product
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent
    / "open_problem_forr4_passive_floor_consolidation_round_3"
    / "searches"
)
sys.path.insert(0, str(ROOT / "searches"))
sys.path.insert(0, str(ROUND3_SEARCHES))

from q64_noncubic_recovered_universal_insertion import (  # noqa: E402
    REMAINING_CLASS_LABELS,
    artifact_text,
    character_product_factor,
    degree_seven_entries,
    diagnostic,
    inserted_coefficients,
    noncubic_coefficient_map,
    noncubic_entries,
    recovered_universal_entries,
    universal_noncubic_entries,
)
from q64_paper_target_gate import (  # noqa: E402
    RESERVE_TARGET,
    THRESHOLD,
    optimize,
)
from q64_remaining_class_gates import partition_remaining  # noqa: E402


def character(left: int, right: int) -> int:
    return -1 if int(left & right).bit_count() % 2 else 1


def injective_character_average(
    order: int,
    labels: tuple[int, ...],
    excluded: int,
) -> float:
    values = tuple(value for value in range(order) if value != excluded)
    total = 0
    count = 0
    for images in permutations(values, len(labels)):
        term = 1
        for label, image in zip(labels, images, strict=True):
            term *= character(label, image)
        total += term
        count += 1
    return total / count


def exact_character_maxima(order: int) -> tuple[float, float, float]:
    nonzero = tuple(range(1, order))
    maxima = []
    for size in (1, 2, 3):
        maximum = 0.0
        for labels in product(nonzero, repeat=size):
            for excluded in range(order):
                maximum = max(
                    maximum,
                    abs(
                        injective_character_average(
                            order, labels, excluded
                        )
                    ),
                )
        maxima.append(maximum)
    return tuple(maxima)


def parity_and_active_groups(
    order: int,
    support: tuple[int, ...],
    axis: int,
) -> tuple[int, int]:
    groups: dict[int, list[int]] = {}
    for coordinate in support:
        labels = divmod(coordinate, order)
        groups.setdefault(labels[axis], []).append(labels[1 - axis])
    odd = 0
    active = 0
    for values in groups.values():
        if len(values) % 2:
            odd += 1
            continue
        xor = 0
        for value in values:
            xor ^= value
        active += xor != 0
    return odd, active


def direct_q4_endpoint_product_maximum() -> float:
    order = 4
    dimension = order * order
    hadamard = np.asarray([[1]], dtype=np.int8)
    while len(hadamard) < order:
        hadamard = np.block(
            [[hadamard, hadamard], [hadamard, -hadamard]]
        )
    left = []
    right = []
    for permutation in permutations(range(order)):
        for signs in product((-1, 1), repeat=order):
            signed_permutation = np.zeros(
                (order, order), dtype=np.int8
            )
            for column, row in enumerate(permutation):
                signed_permutation[row, column] = signs[column]
            left.append((hadamard @ signed_permutation).reshape(-1))
            right.append((signed_permutation @ hadamard).reshape(-1))
    left_array = np.asarray(left, dtype=np.int8)
    right_array = np.asarray(right, dtype=np.int8)
    supports = tuple(combinations(range(dimension), 7))
    support_array = np.asarray(supports, dtype=np.int16)
    left_features = np.prod(
        left_array[:, support_array], axis=2, dtype=np.int8
    )
    right_features = np.prod(
        right_array[:, support_array], axis=2, dtype=np.int8
    )
    normalization = len(left_array)
    endpoint_forward = (
        left_array.T.astype(float) @ right_features / normalization
    )
    endpoint_reverse = (
        left_features.T.astype(float) @ right_array / normalization
    )
    maximum = 0.0
    for index, support in enumerate(supports):
        row_record, _ = parity_and_active_groups(order, support, 0)
        column_record, _ = parity_and_active_groups(order, support, 1)
        if row_record != 1 or column_record != 1:
            continue
        maximum = max(
            maximum,
            float(np.max(np.abs(endpoint_forward[:, index])))
            * float(np.max(np.abs(endpoint_reverse[index, :]))),
        )
    return maximum


def main() -> None:
    result = diagnostic()
    noncubic = noncubic_entries()
    degree_seven = degree_seven_entries()
    universal_noncubic = universal_noncubic_entries()
    recovered = recovered_universal_entries()
    if len(noncubic) != 140 or set(noncubic) != set(
        partition_remaining()["noncubic_profile"]
    ):
        raise AssertionError("noncubic inventory")
    if len(degree_seven) != 16 or len(universal_noncubic) != 124:
        raise AssertionError("noncubic theorem partition")
    if len(recovered) != 96:
        raise AssertionError("recovered universal inventory")
    if set(noncubic_coefficient_map()) != set(noncubic):
        raise AssertionError("noncubic coefficient coverage")

    for order in (4, 8):
        observed = exact_character_maxima(order)
        expected = (
            1 / (order - 1),
            1 / (order - 1),
            3 / ((order - 1) * (order - 3)),
        )
        if not np.allclose(observed, expected, rtol=0, atol=3e-14):
            raise AssertionError(("injective character maxima", order, observed))

    graph_pairs: set[tuple[int, int]] = set()
    for support in combinations(range(16), 7):
        row_record, row_active = parity_and_active_groups(4, support, 0)
        column_record, column_active = parity_and_active_groups(
            4, support, 1
        )
        if row_record != 1 or column_record != 1:
            continue
        graph_pairs.add((row_active, column_active))
        if row_active == 0 and column_active != 3:
            raise AssertionError(("zero-active row lemma", support))
        if column_active == 0 and row_active != 3:
            raise AssertionError(("zero-active column lemma", support))
    expected_pairs = {
        (0, 3),
        (1, 3),
        (2, 2),
        (2, 3),
        (3, 0),
        (3, 1),
        (3, 2),
        (3, 3),
    }
    if graph_pairs != expected_pairs:
        raise AssertionError(("q4 active-group pairs", graph_pairs))
    direct_product = direct_q4_endpoint_product_maximum()
    if not np.isclose(
        direct_product,
        character_product_factor(4) / 4**2,
        rtol=0,
        atol=3e-14,
    ):
        raise AssertionError(("q4 endpoint-product maximum", direct_product))

    discrete = (
        result.degree_seven_entries,
        result.degree_seven_orbits,
        result.universal_noncubic_entries,
        result.noncubic_entries,
        result.recovered_universal_entries,
        result.previous_proved_entries,
        result.post_noncubic_proved_entries,
        result.total_proved_entries,
        result.remaining_open_entries,
        result.remaining_class_counts,
    )
    if discrete != (16, 4, 124, 140, 96, 428, 568, 664, 224, (176, 48)):
        raise AssertionError(("noncubic/recovered discrete", discrete))
    observed = (
        result.injective_character_product_factor,
        result.degree_seven_minimum_coefficient,
        result.degree_seven_maximum_coefficient,
        result.recovered_universal_coefficient,
        result.previous_routing.total,
        result.post_noncubic_routing.total,
        result.post_noncubic_routing.margin_to_one_third,
        result.final_routing.total,
        result.final_routing.margin_to_one_third,
        result.noncubic_routing_improvement,
        result.recovered_universal_routing_cost,
        result.net_routing_change,
        result.reserve_after_declared_allowance,
        result.adaptive_multiplier_cap_retaining_allowance,
    )
    expected = (
        0.00078064012490242,
        0.01708990969034196,
        0.03823058831532936,
        1.0,
        0.3238115631713356,
        0.2869020767941879,
        0.04643125653914543,
        0.3289382301229411,
        0.004395103210392215,
        0.0369094863771477,
        0.0420361533287532,
        0.0051266669516055186,
        0.003395103210392215,
        1.0103214004924976,
    )
    if not np.allclose(observed, expected, rtol=3e-9, atol=3e-12):
        raise AssertionError(("noncubic/recovered numeric", observed))

    partition = partition_remaining()
    base = inserted_coefficients()
    for label, count, frozen, gate in zip(
        result.remaining_class_labels,
        result.remaining_class_counts,
        result.remaining_class_frozen_targets,
        result.remaining_class_reserve_gates,
        strict=True,
    ):
        if label not in REMAINING_CLASS_LABELS:
            raise AssertionError(("remaining class label", label))
        entries = partition[label]
        if len(entries) != count:
            raise AssertionError(("remaining class count", label))
        if not all(base[entry] == frozen for entry in entries):
            raise AssertionError(("remaining frozen target", label))
        trial = dict(base)
        for entry in entries:
            trial[entry] = gate
        gate_total = optimize(mapped_coefficients=trial).total
        if not np.isclose(
            gate_total,
            THRESHOLD - RESERVE_TARGET,
            rtol=0,
            atol=4e-10,
        ):
            raise AssertionError(("remaining reserve gate", label, gate_total))

    committed = (
        ROOT
        / "artifacts"
        / "q64_noncubic_recovered_universal_insertion.json"
    ).read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale noncubic/recovered insertion artifact")
    print(
        "q64 noncubic/recovered universal insertion passed: "
        f"noncubic={result.noncubic_entries},"
        f"degree7_range={result.degree_seven_minimum_coefficient:.12g}/"
        f"{result.degree_seven_maximum_coefficient:.12g},"
        f"recovered={result.recovered_universal_entries},"
        f"proved_entries={result.total_proved_entries},"
        f"total={result.final_routing.total:.12g},"
        f"margin={result.final_routing.margin_to_one_third:.12g},"
        f"remaining_open={result.remaining_open_entries},"
        "status=mixed_16_proved_220_quarantined_cumulative_values_withdrawn"
    )


if __name__ == "__main__":
    main()
