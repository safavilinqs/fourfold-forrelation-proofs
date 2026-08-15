#!/usr/bin/env python3
"""Regression for the all-projective repair of the terminal three-path witness."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from terminal_three_path_projective_repair import (  # noqa: E402
    PATH_VERTICES,
    flattening_exponent,
    grouped_injective_upper_exponent,
    maximum_occupancy,
    path_partitions,
    three_path_projective_certificate,
)


def sylvester(dimension: int) -> np.ndarray:
    """Return the normalized Sylvester matrix of the requested dimension."""

    hadamard = np.array([[1.0]])
    while hadamard.shape[0] < dimension:
        hadamard = np.block([[hadamard, hadamard], [hadamard, -hadamard]])
    return hadamard / np.sqrt(dimension)


def path_tensor(dimension: int) -> np.ndarray:
    """Return H_ab H_bc H_cd as a four-index tensor."""

    hadamard = sylvester(dimension)
    return np.einsum("ab,bc,cd->abcd", hadamard, hadamard, hadamard)


def direct_flattening_checks() -> None:
    """Compare every path cut with the exact cut-rank formula."""

    for dimension in (2, 4):
        tensor = path_tensor(dimension)
        for size in range(1, len(PATH_VERTICES)):
            for selected_tuple in combinations(PATH_VERTICES, size):
                selected = frozenset(selected_tuple)
                complement = tuple(
                    vertex for vertex in PATH_VERTICES if vertex not in selected
                )
                axes = selected_tuple + complement
                rows = dimension ** len(selected_tuple)
                columns = dimension ** len(complement)
                matrix = np.transpose(tensor, axes).reshape(rows, columns)
                observed = float(np.linalg.norm(matrix, ord=2))
                expected = dimension ** float(flattening_exponent(selected))
                if not np.isclose(observed, expected, atol=3e-12):
                    raise AssertionError(
                        (
                            "path flattening formula",
                            dimension,
                            selected,
                            observed,
                            expected,
                        )
                    )


def main() -> None:
    partitions = path_partitions()
    if len(partitions) != 15:
        raise AssertionError(("Bell number B4", len(partitions)))
    occupancy_counts = {
        occupancy: sum(
            maximum_occupancy(partition) == occupancy for partition in partitions
        )
        for occupancy in range(1, 5)
    }
    if occupancy_counts != {1: 1, 2: 9, 3: 4, 4: 1}:
        raise AssertionError(("partition occupancy classes", occupancy_counts))

    singleton = ((0,), (1,), (2,), (3,))
    if grouped_injective_upper_exponent(singleton) != Fraction(-1, 2):
        raise AssertionError("all-singleton path should retain N^{-1/2}")
    adjacent_pair = ((0, 1), (2,), (3,))
    if grouped_injective_upper_exponent(adjacent_pair) != 0:
        raise AssertionError("the explicit adjacent pair should cost at most one")
    for partition in partitions:
        if maximum_occupancy(partition) == 2:
            exponent = grouped_injective_upper_exponent(partition)
            if exponent > 0:
                raise AssertionError(
                    ("two-mark path amplification", partition, exponent)
                )

    certificate = three_path_projective_certificate()
    if certificate.strong_worst_exponent != 0:
        raise AssertionError(("strong path exponent", certificate))
    if certificate.weak_exponent != Fraction(-1, 2):
        raise AssertionError(("weak path exponent", certificate))
    if certificate.combined_exponent != -1:
        raise AssertionError(("three-path exponent", certificate))
    if certificate.extra_gain_beyond_accepted != Fraction(-1, 2):
        raise AssertionError(("gain beyond global dichotomy", certificate))
    if certificate.target_slack != Fraction(-1, 4):
        raise AssertionError(("level-twelve target slack", certificate))
    if certificate.distinctness_pairs != 12:
        raise AssertionError(("same-layer mask pairs", certificate))

    direct_flattening_checks()

    dimension = 1024
    print(
        "terminal three-path projective repair passed: "
        f"partitions={certificate.partition_count},"
        f"occupancy_classes={occupancy_counts},"
        f"strong_worst={certificate.strong_worst_exponent},"
        f"weak={certificate.weak_exponent},"
        f"combined={certificate.combined_exponent},"
        f"accepted={certificate.accepted_exponent},"
        f"target={certificate.level_twelve_target_exponent},"
        f"N1024_bound={dimension ** float(certificate.combined_exponent):.12g},"
        f"N1024_target={dimension ** float(certificate.level_twelve_target_exponent):.12g}"
    )


if __name__ == "__main__":
    main()
