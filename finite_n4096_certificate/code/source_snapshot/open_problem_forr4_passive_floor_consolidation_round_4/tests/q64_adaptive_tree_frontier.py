#!/usr/bin/env python3
"""Independent algebra and witness checks for the adaptive tree frontier."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from json import loads
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_adaptive_tree_frontier import artifact_text, certificate  # noqa: E402


K4 = (
    (1, 1, 1, 1),
    (1, -1, 1, -1),
    (1, 1, -1, -1),
    (1, -1, -1, 1),
)


def exact_two_copy_witness() -> None:
    """Reconstruct the old square-function witness inside the new frontier."""

    def probability(outcome: int, signs: tuple[int, ...]) -> Fraction:
        amplitude = sum(K4[outcome][block] * signs[block] for block in range(4))
        return Fraction(amplitude * amplitude, 16)

    signed_output = tuple(
        tuple(
            sum(
                (
                    Fraction(signs[0] * signs[1] * signs[2] * signs[3], 16)
                    * probability(first, signs)
                    * probability(second, signs)
                )
                for signs in product((-1, 1), repeat=4)
            )
            for second in range(4)
        )
        for first in range(4)
    )
    expected = tuple(
        tuple(Fraction(3 if row == column else -1, 32) for column in range(4))
        for row in range(4)
    )
    if signed_output != expected:
        raise AssertionError((signed_output, expected))
    transcript_mass = sum(abs(value) for row in signed_output for value in row)
    if transcript_mass != Fraction(3, 4):
        raise AssertionError(transcript_mass)

    effects = tuple(
        tuple(
            tuple(Fraction(K4[y][i] * K4[y][j], 4) for i in range(4))
            for j in range(4)
        )
        for y in range(4)
    )
    observables = []
    for y in range(4):
        observables.append(
            tuple(
                tuple(
                    sum(
                        (Fraction(1 if y == z else -1) * effects[z][ell][k]
                         for z in range(4)),
                        Fraction(0),
                    )
                    for k in range(4)
                )
                for ell in range(4)
            )
        )

    # D is the complete signed two-batch strategy kernel, with both uniform
    # probe amplitudes included. Its exact feature masses are 1/16 per
    # temporal row and column history.
    strategy = {}
    for i, k, j, ell in product(range(4), repeat=4):
        strategy[i, k, j, ell] = Fraction(1, 16) * sum(
            (effects[y][j][i] * observables[y][ell][k] for y in range(4)),
            Fraction(0),
        )

    # Build the direct-sum features explicitly. For a rank-one Hadamard
    # effect, sqrt(E_y)=E_y, and H_y=2E_y-I is an orthogonal reflection.
    row_features = {}
    column_features = {}
    for i, k in product(range(4), repeat=2):
        row = []
        column = []
        for y, root_coordinate, child_coordinate in product(range(4), repeat=3):
            root_atom = Fraction(K4[y][root_coordinate] * K4[y][i], 8)
            row.append(root_atom * Fraction(1, 2) * observables[y][child_coordinate][k])
            column.append(
                root_atom * Fraction(1, 2)
                * Fraction(1 if child_coordinate == k else 0)
            )
        row_features[i, k] = tuple(row)
        column_features[i, k] = tuple(column)

    for i, k, j, ell in product(range(4), repeat=4):
        inner = sum(
            (
                row_features[i, k][coordinate]
                * column_features[j, ell][coordinate]
                for coordinate in range(64)
            ),
            Fraction(0),
        )
        if inner != strategy[i, k, j, ell]:
            raise AssertionError((i, k, j, ell, inner, strategy[i, k, j, ell]))
    for features in (row_features, column_features):
        masses = {
            history: sum((value * value for value in vector), Fraction(0))
            for history, vector in features.items()
        }
        if set(masses.values()) != {Fraction(1, 16)} or sum(masses.values()) != 1:
            raise AssertionError(masses)

    def signed_moment(i: int, k: int, j: int, ell: int) -> Fraction:
        return sum(
            (
                Fraction(signs[0] * signs[1] * signs[2] * signs[3], 16)
                * signs[i] * signs[k] * signs[j] * signs[ell]
                for signs in product((-1, 1), repeat=4)
            ),
            Fraction(0),
        )

    reconstructed = sum(
        (
            strategy[i, k, j, ell] * signed_moment(i, k, j, ell)
            for i, k, j, ell in product(range(4), repeat=4)
        ),
        Fraction(0),
    )
    if reconstructed != transcript_mass:
        raise AssertionError((reconstructed, transcript_mass))

    # The rejected proxy used joint square mass 3/32 and therefore demanded
    # sqrt(6) at N=1. The normalized strategy factorization above has mass one.
    if transcript_mass * transcript_mass / Fraction(3, 32) != 6:
        raise AssertionError("old square-function witness changed")


def random_recursive_factorization() -> None:
    """Check the direct-sum induction with arbitrary complex POVMs."""

    generator = np.random.default_rng(2026071819)
    for _ in range(80):
        root_labels, idler, outcomes = 3, 2, 4
        root_dimension = root_labels * idler
        raw_amplitudes = (
            generator.normal(size=(root_labels, idler))
            + 1j * generator.normal(size=(root_labels, idler))
        )
        raw_amplitudes /= np.linalg.norm(raw_amplitudes)
        psi = np.zeros((root_labels, root_dimension), dtype=complex)
        for label in range(root_labels):
            psi[label, label * idler : (label + 1) * idler] = raw_amplitudes[label]

        raw_kraus = (
            generator.normal(size=(outcomes, root_dimension, root_dimension))
            + 1j * generator.normal(size=(outcomes, root_dimension, root_dimension))
        )
        gram = sum((matrix.conj().T @ matrix for matrix in raw_kraus))
        eigenvalues, eigenvectors = np.linalg.eigh(gram)
        inverse_root = (eigenvectors * (1 / np.sqrt(eigenvalues))) @ eigenvectors.conj().T
        kraus = np.einsum("yab,bc->yac", raw_kraus, inverse_root)
        completeness = sum((matrix.conj().T @ matrix for matrix in kraus))
        if not np.allclose(completeness, np.eye(root_dimension), atol=3e-12):
            raise AssertionError("POVM normalization")
        root_features = np.einsum("yab,ib->yia", kraus, psi)

        child_rows, child_columns, child_dimension = 3, 4, 5
        child_left = (
            generator.normal(size=(outcomes, child_rows, child_dimension))
            + 1j * generator.normal(size=(outcomes, child_rows, child_dimension))
        )
        child_right = (
            generator.normal(size=(outcomes, child_columns, child_dimension))
            + 1j * generator.normal(size=(outcomes, child_columns, child_dimension))
        )
        for y in range(outcomes):
            child_left[y] /= np.linalg.norm(child_left[y])
            child_right[y] /= np.linalg.norm(child_right[y])

        row_features = np.empty(
            (root_labels * child_rows, outcomes * root_dimension * child_dimension),
            dtype=complex,
        )
        column_features = np.empty(
            (root_labels * child_columns, outcomes * root_dimension * child_dimension),
            dtype=complex,
        )
        for i, child in product(range(root_labels), range(child_rows)):
            row_features[i * child_rows + child] = np.concatenate(
                [np.kron(root_features[y, i], child_left[y, child]) for y in range(outcomes)]
            )
        for j, child in product(range(root_labels), range(child_columns)):
            column_features[j * child_columns + child] = np.concatenate(
                [np.kron(root_features[y, j], child_right[y, child]) for y in range(outcomes)]
            )

        direct = np.zeros((root_labels * child_rows, root_labels * child_columns), dtype=complex)
        for i, child_i, j, child_j in product(
            range(root_labels), range(child_rows), range(root_labels), range(child_columns)
        ):
            direct[i * child_rows + child_i, j * child_columns + child_j] = sum(
                (
                    np.vdot(root_features[y, j], root_features[y, i])
                    * np.vdot(child_right[y, child_j], child_left[y, child_i])
                    for y in range(outcomes)
                )
            )
        factored = row_features @ column_features.conj().T
        if not np.allclose(direct, factored, atol=5e-12):
            raise AssertionError("direct-sum factorization")
        if np.linalg.norm(row_features) ** 2 > 1 + 3e-12:
            raise AssertionError("row feature mass")
        if np.linalg.norm(column_features) ** 2 > 1 + 3e-12:
            raise AssertionError("column feature mass")


def rare_outcome_factorization() -> None:
    """Exercise the induction without dividing by a rare outcome mass."""

    epsilon = 1e-14
    amplitudes = np.sqrt(np.array([0.6, 0.4]))
    effects = (
        np.diag([epsilon, epsilon * epsilon]),
        np.diag([1 - epsilon, 1 - epsilon * epsilon]),
    )
    child_left = (
        np.array([[1.0, 0.0], [0.0, 1.0]]) / np.sqrt(2),
        np.array([[1.0, 1.0], [1.0, -1.0]]) / 2,
    )
    child_right = (
        np.array([[1.0, 1.0], [1.0, -1.0]]) / 2,
        np.array([[1.0, 0.0], [0.0, 1.0]]) / np.sqrt(2),
    )
    root_features = np.zeros((2, 2, 2))
    for y, effect in enumerate(effects):
        root_features[y] = np.diag(np.sqrt(np.diag(effect))) * amplitudes[:, None]

    rows = []
    columns = []
    for root, child in product(range(2), repeat=2):
        rows.append(
            np.concatenate(
                [np.kron(root_features[y, root], child_left[y][child]) for y in range(2)]
            )
        )
        columns.append(
            np.concatenate(
                [np.kron(root_features[y, root], child_right[y][child]) for y in range(2)]
            )
        )
    row_features = np.asarray(rows)
    column_features = np.asarray(columns)
    if np.linalg.norm(row_features) ** 2 > 1 + 2e-14:
        raise AssertionError("rare row mass")
    if np.linalg.norm(column_features) ** 2 > 1 + 2e-14:
        raise AssertionError("rare column mass")
    rare_mass = sum(
        amplitudes[root] ** 2 * effects[0][root, root]
        for root in range(2)
    )
    if not (0 < rare_mass < 7e-15):
        raise AssertionError(("rare outcome", rare_mass))


def dose_partition_pullbacks() -> None:
    """Check aggregate characters for all lead two-batch dose splits."""

    generator = np.random.default_rng(2026071821)
    modes = 12
    for root_dose, child_dose in ((1, 5), (2, 4), (3, 3)):
        for _ in range(120):
            root = np.bincount(
                generator.integers(0, modes, size=root_dose), minlength=modes
            )
            child = np.bincount(
                generator.integers(0, modes, size=child_dose), minlength=modes
            )
            signs = generator.choice((-1, 1), size=modes)
            root_character = int(np.prod(signs ** root))
            child_character = int(np.prod(signs ** child))
            aggregate_character = int(np.prod(signs ** (root + child)))
            if root_character * child_character != aggregate_character:
                raise AssertionError((root_dose, child_dose, root, child, signs))
            if int((root + child).sum()) != 6:
                raise AssertionError("branch dose")


def pullback_and_perron_checks() -> None:
    """Stress duplicate histories and independent occupation laws."""

    generator = np.random.default_rng(2026071820)
    for _ in range(100):
        base_rows, base_columns, base_rank = 4, 5, 3
        left = generator.normal(size=(base_rows, base_rank))
        right = generator.normal(size=(base_columns, base_rank))
        left /= np.maximum(np.linalg.norm(left, axis=1, keepdims=True), 1.0)
        right /= np.maximum(np.linalg.norm(right, axis=1, keepdims=True), 1.0)
        histories_left, histories_right, temporal_rank = 8, 9, 4
        row_map = generator.integers(0, base_rows, size=histories_left)
        column_map = generator.integers(0, base_columns, size=histories_right)
        temporal_left = generator.normal(size=(histories_left, temporal_rank))
        temporal_right = generator.normal(size=(histories_right, temporal_rank))
        temporal_left /= np.linalg.norm(temporal_left)
        temporal_right /= np.linalg.norm(temporal_right)
        matrix = (left[row_map] @ right[column_map].T) * (
            temporal_left @ temporal_right.T
        )
        nuclear = float(np.linalg.svd(matrix, compute_uv=False).sum())
        factor_bound = float(
            np.max(np.linalg.norm(left, axis=1))
            * np.max(np.linalg.norm(right, axis=1))
        )
        if nuclear > factor_bound + 2e-12:
            raise AssertionError(("pullback factorization", nuclear, factor_bound))

        size = 7
        raw = generator.random((size, size))
        perron_matrix = (raw + raw.T) / 2
        spectral = float(np.linalg.eigvalsh(perron_matrix)[-1])
        first = generator.random(size)
        second = generator.random(size)
        first /= np.linalg.norm(first)
        second /= np.linalg.norm(second)
        bilinear = float(first @ perron_matrix @ second)
        if bilinear > spectral + 2e-12:
            raise AssertionError(("independent Perron laws", bilinear, spectral))


def main() -> None:
    result = certificate()
    if result.verdict != "CERTIFIED" or not result.passes_reserve_gate:
        raise AssertionError(result)
    if (result.frontier_multiplier_numerator, result.frontier_multiplier_denominator) != (1, 1):
        raise AssertionError("frontier multiplier")
    exact_two_copy_witness()
    random_recursive_factorization()
    rare_outcome_factorization()
    dose_partition_pullbacks()
    pullback_and_perron_checks()

    committed_path = ROOT / "artifacts" / "q64_adaptive_tree_frontier.json"
    committed = committed_path.read_text(encoding="utf-8")
    if committed != artifact_text(result):
        raise AssertionError("stale adaptive-tree artifact")
    payload = loads(committed)
    if payload["result"]["verdict"] != "CERTIFIED":
        raise AssertionError("artifact verdict")
    print(
        "q64 adaptive tree frontier passed: "
        "two_copy=exact,recursive_factorizations=80,"
        "rare_outcome_below=7e-15,dose_partitions=1+5/2+4/3+3,"
        "pullback_perron_checks=100,multiplier=1,"
        f"total={result.adaptive_total_upper},verdict=certified"
    )


if __name__ == "__main__":
    main()
