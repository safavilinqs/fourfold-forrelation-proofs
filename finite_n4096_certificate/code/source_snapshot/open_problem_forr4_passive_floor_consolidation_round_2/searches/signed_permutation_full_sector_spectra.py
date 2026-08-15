#!/usr/bin/env python3
"""Full low-degree signed-permutation link sectors, including decorations."""

from __future__ import annotations

from itertools import combinations, permutations, product
from math import comb

import numpy as np


Q = 4
N = Q * Q


def sylvester_sign(order: int) -> np.ndarray:
    result = np.array([[1]], dtype=np.int8)
    while result.shape[0] < order:
        result = np.block([[result, result], [result, -result]])
    return result


def all_pairs() -> tuple[np.ndarray, np.ndarray]:
    k = sylvester_sign(Q)
    left = []
    right = []
    for permutation in permutations(range(Q)):
        for signs in product((-1, 1), repeat=Q):
            p = np.zeros((Q, Q), dtype=np.int8)
            for column, row in enumerate(permutation):
                p[row, column] = signs[column]
            left.append((k @ p).reshape(-1))
            right.append((p @ k).reshape(-1))
    return np.array(left, dtype=np.int8), np.array(right, dtype=np.int8)


def supports(degree: int) -> list[tuple[int, ...]]:
    return list(combinations(range(N), degree))


def features(values: np.ndarray, sector: list[tuple[int, ...]]) -> np.ndarray:
    result = np.empty((len(values), len(sector)), dtype=np.int8)
    for column, support in enumerate(sector):
        result[:, column] = np.prod(values[:, support], axis=1)
    return result


def record_size(support: tuple[int, ...], side: str) -> int:
    counts = np.zeros(Q, dtype=np.int8)
    for coordinate in support:
        row, column = divmod(coordinate, Q)
        label = column if side == "left" else row
        counts[label] ^= 1
    return int(counts.sum())


def one_label_type(support: tuple[int, ...], side: str) -> bool:
    labels = []
    for coordinate in support:
        row, column = divmod(coordinate, Q)
        labels.append(column if side == "left" else row)
    return len(set(labels)) == 1


def operator_norm_from_features(left: np.ndarray, right: np.ndarray) -> float:
    count = len(left)
    left_gram = left.astype(float) @ left.T / count
    right_gram = right.astype(float) @ right.T / count
    left_values, left_vectors = np.linalg.eigh(left_gram)
    right_values, right_vectors = np.linalg.eigh(right_gram)
    left_keep = left_values > 1e-10
    right_keep = right_values > 1e-10
    core = (
        np.sqrt(left_values[left_keep])[:, None]
        * (left_vectors[:, left_keep].T @ right_vectors[:, right_keep])
        * np.sqrt(right_values[right_keep])[None, :]
    )
    return float(np.linalg.svd(core, compute_uv=False)[0])


def main() -> None:
    left_values, right_values = all_pairs()
    support_families = {degree: supports(degree) for degree in (1, 3)}
    left_features = {
        degree: features(left_values, support_families[degree])
        for degree in support_families
    }
    right_features = {
        degree: features(right_values, support_families[degree])
        for degree in support_families
    }

    rows = []
    for left_degree, right_degree in product((1, 3), repeat=2):
        left_sector = support_families[left_degree]
        right_sector = support_families[right_degree]
        left_records = np.array(
            [record_size(support, "left") for support in left_sector]
        )
        right_records = np.array(
            [record_size(support, "right") for support in right_sector]
        )
        for record in sorted(set(left_records) & set(right_records)):
            left_mask = left_records == record
            right_mask = right_records == record
            lf = left_features[left_degree][:, left_mask]
            rf = right_features[right_degree][:, right_mask]
            moment = lf.T.astype(float) @ rf / len(left_values)
            op = operator_norm_from_features(lf, rf)
            maximum = float(np.max(np.abs(moment)))
            row_energy = float(np.max(np.sum(moment * moment, axis=1)))
            column_energy = float(np.max(np.sum(moment * moment, axis=0)))
            rows.append(
                f"degrees={left_degree},{right_degree},record={record},"
                f"shape={moment.shape[0]}x{moment.shape[1]},"
                f"max_entry={maximum:.12g},op={op:.12g},"
                f"max_row_energy={row_energy:.12g},"
                f"max_column_energy={column_energy:.12g}"
            )

            expected_entry = 1 / comb(Q, record)
            if maximum > expected_entry * (1 + 2e-12):
                raise AssertionError(
                    ("matching coherence", left_degree, right_degree, record)
                )

    m13 = (
        left_features[1].T.astype(float) @ right_features[3]
        / len(left_values)
    )
    m31 = (
        left_features[3].T.astype(float) @ right_features[1]
        / len(left_values)
    )
    m11 = (
        left_features[1].T.astype(float) @ right_features[1]
        / len(left_values)
    )
    expected_full_energy = (Q * Q + 2) / 6
    active_m13 = m13[:, np.any(np.abs(m13) > 2e-12, axis=0)]
    full_gram = active_m13 @ active_m13.T
    if not np.allclose(
        full_gram,
        expected_full_energy * np.eye(N),
        atol=2e-12,
    ):
        raise AssertionError(("full cubic decorated Gram", full_gram))
    profile = (
        m13[:, :, None, None]
        * m31[None, :, :, None]
        * m11[None, None, :, :]
    )
    alternating = np.transpose(profile, (0, 2, 1, 3)).reshape(
        m13.shape[0] * m31.shape[1],
        m13.shape[1] * m11.shape[1],
    )
    uniform_nuclear = float(
        np.linalg.svd(alternating, compute_uv=False).sum()
        / np.sqrt(alternating.shape[0] * alternating.shape[1])
    )
    rows.append(
        "profile=1,3,1,1,cut=13|24,"
        f"shape={alternating.shape[0]}x{alternating.shape[1]},"
        f"uniform_normalized_nuclear={uniform_nuclear:.12g}"
    )
    middle_supports = support_families[3]
    compatible_middle = sum(
        record_size(support, "left") == 1
        and record_size(support, "right") == 1
        for support in middle_supports
    )
    expected_middle = Q * Q * (Q - 1) * (Q - 1)
    if compatible_middle != expected_middle:
        raise AssertionError(
            ("aligned middle support count", compatible_middle, expected_middle)
        )
    compatible_mask = np.array(
        [
            record_size(support, "left") == 1
            and record_size(support, "right") == 1
            for support in middle_supports
        ]
    )
    restricted13 = m13[:, compatible_mask]
    restricted31 = m31[compatible_mask, :]
    for name, matrix in (
        ("M13", restricted13),
        ("M31", restricted31),
    ):
        restricted_op = float(np.linalg.svd(matrix, compute_uv=False)[0])
        restricted_entry = float(np.max(np.abs(matrix)))
        if not np.isclose(restricted_op, 1.0, atol=2e-12):
            raise AssertionError(("compatible L operator", name, restricted_op))
        if not np.isclose(
            restricted_entry, 1 / (Q * (Q - 1)), atol=2e-12
        ):
            raise AssertionError(
                ("compatible L coherence", name, restricted_entry)
            )
    rows.append(
        f"profile=1,3,1,1,compatible_middle_supports={compatible_middle},"
        "restricted_outer_op=1,"
        f"restricted_coherence={1 / (Q * (Q - 1)):.12g}"
    )

    endpoint_profile = (
        m31[:, :, None, None]
        * m11[None, :, :, None]
        * m11[None, None, :, :]
    )
    endpoint_alternating = np.transpose(
        endpoint_profile, (0, 2, 1, 3)
    ).reshape(
        m31.shape[0] * m11.shape[0],
        m31.shape[1] * m11.shape[1],
    )
    endpoint_uniform_nuclear = float(
        np.linalg.svd(endpoint_alternating, compute_uv=False).sum()
        / np.sqrt(
            endpoint_alternating.shape[0] * endpoint_alternating.shape[1]
        )
    )
    rows.append(
        "profile=3,1,1,1,cut=13|24,"
        f"shape={endpoint_alternating.shape[0]}x"
        f"{endpoint_alternating.shape[1]},"
        "uniform_normalized_nuclear="
        f"{endpoint_uniform_nuclear:.12g}"
    )

    endpoint_type_mask = np.array(
        [one_label_type(support, "left") for support in support_families[3]]
    )
    type_one_m31 = m31[endpoint_type_mask, :]
    endpoint_type_profile = (
        type_one_m31[:, :, None, None]
        * m11[None, :, :, None]
        * m11[None, None, :, :]
    )
    endpoint_type_alternating = np.transpose(
        endpoint_type_profile, (0, 2, 1, 3)
    ).reshape(
        type_one_m31.shape[0] * m11.shape[0],
        type_one_m31.shape[1] * m11.shape[1],
    )
    endpoint_type_uniform = float(
        np.linalg.svd(endpoint_type_alternating, compute_uv=False).sum()
        / np.sqrt(
            endpoint_type_alternating.shape[0]
            * endpoint_type_alternating.shape[1]
        )
    )
    if not np.isclose(endpoint_type_uniform, 1 / Q, atol=2e-12):
        raise AssertionError(
            ("endpoint cubic normalized nuclear", endpoint_type_uniform)
        )
    rows.append(
        "profile=3,1,1,1,endpoint_one_label_orbit="
        f"{endpoint_type_alternating.shape[0]}x"
        f"{endpoint_type_alternating.shape[1]},"
        f"uniform_normalized_nuclear={endpoint_type_uniform:.12g}"
    )

    print("signed-permutation full sectors:\n" + "\n".join(rows))


if __name__ == "__main__":
    main()
