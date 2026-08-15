#!/usr/bin/env python3
"""Attack Case I after enforcing distinct Fourier coordinates.

Two disjoint four-layer Hadamard chains are vertical tensor products before
the Fourier-support constraint is imposed.  Distinctness within each block
multiplies that product by a mask coupling the two components.  This script
computes physical-party flattening upper bounds and higher-order power-method
lower bounds for every relative placement of the second chain.

The exact witness falsifies literal masked multiplicativity.  The numerical
search separately attacks the weaker theorem-level 1/sqrt(N) target; a
passing power search is not a proof of that target.
"""

from __future__ import annotations

from itertools import combinations, permutations

import numpy as np


SEED = 2026071410
TOL = 2e-10


def sylvester(n: int) -> np.ndarray:
    h = np.array([[1.0]])
    while h.shape[0] < n:
        h = np.block([[h, h], [h, -h]])
    return h / np.sqrt(n)


def masked_two_chain_tensor(n: int, placement: tuple[int, ...]) -> np.ndarray:
    """Return the order-four tensor over the four physical entry parties.

    Party e stores (a_e, b_l), where the first-chain layer-e vertex is in
    party e and placement[l] = e for the second-chain layer-l vertex.
    """

    h = sylvester(n)
    party_dim = n * n
    indices = np.arange(party_dim)
    first = indices // n
    second = indices % n

    def axis_factor(left_axis: int, right_axis: int, matrix: np.ndarray) -> np.ndarray:
        """Broadcast matrix[left-coordinate,right-coordinate] on four axes."""

        if left_axis == right_axis:
            shape = [1] * 4
            shape[left_axis] = party_dim
            return np.diag(matrix).reshape(shape)
        low, high = sorted((left_axis, right_axis))
        value = matrix if left_axis < right_axis else matrix.T
        shape = [1] * 4
        shape[low] = party_dim
        shape[high] = party_dim
        return value.reshape(shape)

    first_h = h[first[:, None], first[None, :]]
    second_h = h[second[:, None], second[None, :]]
    tensor: np.ndarray | float = 1.0
    for layer in range(3):
        tensor = tensor * axis_factor(layer, layer + 1, first_h)
        tensor = tensor * axis_factor(placement[layer], placement[layer + 1], second_h)
    for layer in range(4):
        equality = first[:, None] != second[None, :]
        tensor = tensor * axis_factor(layer, placement[layer], equality)
    return np.asarray(tensor, dtype=float)


def flattening_upper_bound(tensor: np.ndarray) -> float:
    """Smallest operator norm over nontrivial physical-party cuts."""

    best = np.inf
    parties = range(tensor.ndim)
    for size in (1, 2):
        for left in combinations(parties, size):
            if size == 2 and 0 not in left:
                continue  # avoid complementary duplicates
            right = tuple(party for party in parties if party not in left)
            order = left + right
            rows = int(np.prod([tensor.shape[party] for party in left]))
            matrix = np.transpose(tensor, order).reshape(rows, -1)
            best = min(best, float(np.linalg.norm(matrix, ord=2)))
    return best


def spectral_lower_bound(
    rng: np.random.Generator, tensor: np.ndarray, restarts: int = 24, sweeps: int = 50
) -> float:
    """Coordinate-ascent lower bound on the complex injective norm."""

    dimension = tensor.shape[0]
    labels = "abcd"
    best = 0.0
    for _ in range(restarts):
        vectors = []
        for _party in range(4):
            vector = rng.normal(size=dimension) + 1j * rng.normal(size=dimension)
            vectors.append(vector / np.linalg.norm(vector))
        previous = 0.0
        for _sweep in range(sweeps):
            for party in range(4):
                other = [index for index in range(4) if index != party]
                expression = labels + "," + ",".join(labels[index] for index in other)
                expression += "->" + labels[party]
                coefficient = np.einsum(
                    expression, tensor, *(vectors[index] for index in other), optimize=True
                )
                norm = np.linalg.norm(coefficient)
                if norm:
                    vectors[party] = np.conj(coefficient) / norm
            value = abs(np.einsum("abcd,a,b,c,d->", tensor, *vectors, optimize=True))
            if abs(value - previous) <= 1e-13 * max(1.0, value):
                break
            previous = value
        best = max(best, float(value))
    return best


def exact_n4_multiplicativity_counterexample() -> float:
    """Return an exact product-vector witness above the unmasked 1/N target."""

    tensor = masked_two_chain_tensor(4, (0, 1, 3, 2))
    first = np.array(
        [
            [0, 1, 1, -1],
            [-1, 0, 1, -1],
            [-1, 1, 0, -1],
            [-1, 1, 1, 0],
        ],
        dtype=float,
    ).ravel() / np.sqrt(12)
    second = np.zeros((4, 4), dtype=float)
    second[0, 3] = -1
    third = np.array(
        [
            [-1, -1, 1, 1],
            [1, -1, -1, 1],
            [-1, -1, -1, -1],
            [1, -1, 1, -1],
        ],
        dtype=float,
    ).ravel() / 4
    fourth = np.array(
        [
            [-1, -1, 1, 1],
            [-1, 1, 1, -1],
            [1, 1, 1, 1],
            [1, -1, 1, -1],
        ],
        dtype=float,
    ).ravel() / 4
    vectors = (first, second.ravel(), third, fourth)
    if not all(np.isclose(np.linalg.norm(vector), 1.0, atol=2e-14) for vector in vectors):
        raise AssertionError("counterexample vectors are not normalized")
    value = abs(np.einsum("abcd,a,b,c,d->", tensor, *vectors, optimize=True))
    expected = 15 / (32 * np.sqrt(3))
    if not np.isclose(value, expected, atol=2e-14):
        raise AssertionError(("exact witness value", value, expected))
    if not value > 1 / 4:
        raise AssertionError(("multiplicativity counterexample disappeared", value))
    return value


def main() -> None:
    rng = np.random.default_rng(SEED)
    exact_witness = exact_n4_multiplicativity_counterexample()
    worst_ratio = 0.0
    worst_placement: tuple[int, ...] = ()
    widest_upper_gap = 0.0
    for n in (2, 4):
        product_target = 1 / n
        theorem_target = 1 / np.sqrt(n)
        for placement in permutations(range(4)):
            tensor = masked_two_chain_tensor(n, placement)
            lower = spectral_lower_bound(rng, tensor)
            upper = flattening_upper_bound(tensor)
            ratio = lower / product_target
            worst_ratio = max(worst_ratio, ratio)
            if ratio == worst_ratio:
                worst_placement = placement
            widest_upper_gap = max(widest_upper_gap, upper / product_target)
            if lower > theorem_target * (1 + TOL):
                raise AssertionError(
                    (
                        "masked tensor violates the theorem-level target",
                        n,
                        placement,
                        lower,
                        theorem_target,
                        upper,
                    )
                )
    print(
        "masked all-singleton multiplicativity defect confirmed: "
        f"exact_witness={exact_witness:.12g}, exact_product_target=0.25, "
        f"worst_search_ratio={worst_ratio:.12g}, placement={worst_placement}, "
        f"widest_flattening_ratio={widest_upper_gap:.12g}"
    )


if __name__ == "__main__":
    main()
