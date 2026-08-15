#!/usr/bin/env python3
"""Regression for the third balanced chain-aware contraction."""

from __future__ import annotations

from itertools import combinations, permutations, product
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from adjacent_cubic_quintic_orbit_witness import (  # noqa: E402
    unnormalized_sylvester,
)
from occupation_compatible_sector_optimization import (  # noqa: E402
    endpoint_quintic_singleton_slice_energies,
    endpoint_singleton_slice_energies,
)
from separated_balanced_endpoint_slice_contraction import (  # noqa: E402
    separated_balanced_coefficient,
    separated_balanced_contraction,
    separated_balanced_orbit_entries,
)


SEED = 2026071503


def exact_q4_endpoint_slices() -> tuple[float, float]:
    """Independently enumerate the M_31 and v_5 fixed-support slices."""

    order = 4
    dimension = order * order
    hadamard_sign = unnormalized_sylvester(order)
    left_values = []
    right_values = []
    for permutation in permutations(range(order)):
        for signs in product((-1, 1), repeat=order):
            signed_permutation = np.zeros((order, order), dtype=np.int8)
            for column, row in enumerate(permutation):
                signed_permutation[row, column] = signs[column]
            left_values.append((hadamard_sign @ signed_permutation).reshape(-1))
            right_values.append((signed_permutation @ hadamard_sign).reshape(-1))
    left = np.asarray(left_values, dtype=np.int8)
    right = np.asarray(right_values, dtype=np.int8)
    triples = tuple(combinations(range(dimension), 3))
    fives = tuple(combinations(range(dimension), 5))
    triple_index = {support: index for index, support in enumerate(triples)}
    five_index = {support: index for index, support in enumerate(fives)}
    left_three = np.asarray(
        [np.prod(left[:, support], axis=1) for support in triples],
        dtype=np.int8,
    ).T
    right_five = np.asarray(
        [np.prod(right[:, support], axis=1) for support in fives],
        dtype=np.int8,
    ).T
    normalization = len(left)
    moment_31 = left_three.T.astype(float) @ right / normalization
    moment_15 = left.T.astype(float) @ right_five / normalization

    maximum_cubic = 0.0
    for singleton in range(dimension):
        for pair in combinations(range(dimension), 2):
            extensions = [
                triple_index[tuple(sorted(pair + (cell,)))]
                for cell in range(dimension)
                if cell not in pair
            ]
            maximum_cubic = max(
                maximum_cubic,
                float(np.sum(moment_31[extensions, singleton] ** 2)),
            )

    normalized_hadamard = np.kron(
        hadamard_sign,
        hadamard_sign,
    ) / order
    maximum_quintic_scalar = 0.0
    for fixed_four in combinations(range(dimension), 4):
        for singleton in range(dimension):
            values = []
            for element in range(dimension):
                if element in fixed_four:
                    continue
                support = tuple(sorted(fixed_four + (element,)))
                moment = moment_15[singleton, five_index[support]]
                xor_support = 0
                for coordinate in support:
                    xor_support ^= coordinate
                hadamard_entry = normalized_hadamard[singleton, xor_support]
                values.append(moment / hadamard_entry)
            maximum_quintic_scalar = max(
                maximum_quintic_scalar,
                float(np.sum(np.asarray(values) ** 2)),
            )
    return maximum_cubic, maximum_quintic_scalar


def cyclic_fourier(dimension: int) -> np.ndarray:
    indices = np.arange(dimension)
    return np.exp(
        2j * np.pi * indices[:, None] * indices[None, :] / dimension
    ) / np.sqrt(dimension)


def direct_abstract_chain_checks() -> float:
    """Stress the separated endpoint factorization on complete tensors."""

    rng = np.random.default_rng(SEED)
    dimension = 6
    pairs = tuple(combinations(range(dimension), 2))
    fours = tuple(combinations(range(dimension), 4))
    fourier = cyclic_fourier(dimension)
    slice_energy = 0.23
    quintic_moment_energy = 0.71
    quintic_scalar_energy = dimension * quintic_moment_energy
    bound = np.sqrt(slice_energy * quintic_moment_energy)
    worst = 0.0

    for _ in range(5):
        cubic = np.zeros((dimension, len(pairs), dimension), dtype=complex)
        for pair_index, pair in enumerate(pairs):
            for singleton in range(dimension):
                values = (
                    rng.normal(size=dimension)
                    + 1j * rng.normal(size=dimension)
                )
                values[list(pair)] = 0
                values *= np.sqrt(slice_energy) / np.linalg.norm(values)
                cubic[:, pair_index, singleton] = values

        quintic = np.zeros((len(fours), dimension), dtype=complex)
        for four_index, four in enumerate(fours):
            values = (
                rng.normal(size=dimension)
                + 1j * rng.normal(size=dimension)
            )
            values[list(four)] = 0
            values *= np.sqrt(quintic_scalar_energy) / np.linalg.norm(values)
            quintic[four_index] = values

        xor_indices = np.asarray(
            [sum(four) % dimension for four in fours],
            dtype=int,
        )
        phases = np.exp(
            2j
            * np.pi
            * rng.random((dimension, dimension, dimension))
        )
        residual = (
            fourier[xor_indices[:, None], np.arange(dimension)[None, :]]
        )[:, None, :, None] * phases[None, :, :, :]

        # Rows are (x,F), columns are (E,b,c,e).  The exact product of the
        # two central Walsh entries contributes 1/sqrt(N) times the residual
        # normalized Fourier Gram symbol.
        tensor = np.einsum(
            "xpb,fe,fbce->xfpbce",
            cubic,
            quintic,
            residual / np.sqrt(dimension),
        ).reshape(
            dimension * len(fours),
            len(pairs) * dimension * dimension * dimension,
        )
        row_law = rng.dirichlet(np.ones(tensor.shape[0]))
        column_law = rng.dirichlet(np.ones(tensor.shape[1]))
        weighted = (
            np.sqrt(row_law)[:, None]
            * tensor
            * np.sqrt(column_law)[None, :]
        )
        nuclear = float(np.linalg.svd(weighted, compute_uv=False).sum())
        if nuclear > bound * (1 + 5e-12):
            raise AssertionError(
                ("abstract separated contraction", nuclear, bound)
            )
        worst = max(worst, nuclear / bound)
    return worst


def main() -> None:
    cubic_slice, quintic_scalar_slice = exact_q4_endpoint_slices()
    expected_cubic = endpoint_singleton_slice_energies(4)[2]
    expected_quintic_moment = endpoint_quintic_singleton_slice_energies(4)[4]
    if not np.isclose(cubic_slice, expected_cubic, atol=2e-14):
        raise AssertionError(("exact M31 fixed-pair slice", cubic_slice))
    if not np.isclose(
        quintic_scalar_slice,
        16 * expected_quintic_moment,
        atol=2e-14,
    ):
        raise AssertionError(
            ("exact v5 fixed-four slice", quintic_scalar_slice)
        )

    orbit = separated_balanced_orbit_entries()
    if len(orbit) != 4:
        raise AssertionError(("separated balanced orbit", orbit))
    coefficient = separated_balanced_coefficient()
    if not np.isclose(coefficient, 0.1737428008469535, atol=2e-15):
        raise AssertionError(("separated balanced coefficient", coefficient))
    result = separated_balanced_contraction()
    if not result.coefficient < result.generic_two_mask_coefficient:
        raise AssertionError(("coefficient does not beat generic masks", result))
    if not result.coefficient > result.provisional_coefficient:
        raise AssertionError(("conservative coefficient relation", result))
    if not result.generic_two_mask_coefficient > 0.225:
        raise AssertionError(("generic obstruction missing", result))
    if not result.threshold_slack > 0.00066:
        raise AssertionError(("three-theorem ledger slack", result))
    if (
        (result.next_unresolved_entries[0])
        != ((3, 1, 1, 5), (1, 1, 1, 2))
    ):
        raise AssertionError(("reranked next orbit", result))
    if not np.isclose(
        result.next_unresolved_contribution,
        0.00150951796806981,
        atol=2e-15,
    ):
        raise AssertionError(("reranked next contribution", result))
    if not np.isclose(
        result.next_admissible_coefficient,
        0.0450405467777823,
        atol=3e-14,
    ):
        raise AssertionError(("next admissible coefficient", result))

    worst = direct_abstract_chain_checks()
    print(
        "separated balanced endpoint-slice contraction passed: "
        f"cubic_slice={result.cubic_pair_slice_energy:.12g},"
        f"quintic_slice={result.quintic_four_slice_energy:.12g},"
        f"coefficient={result.coefficient:.12g},"
        f"ledger_total={result.optimized_total:.12g},"
        f"threshold_slack={result.threshold_slack:.12g},"
        f"next_contribution={result.next_unresolved_contribution:.12g},"
        f"next_admissible={result.next_admissible_coefficient:.12g},"
        f"worst_random_ratio={worst:.12g}"
    )


if __name__ == "__main__":
    main()
