#!/usr/bin/env python3
"""Fixed-seed adaptive passive search for the exact q=2 permutation plant.

The root is a binary rank-one measurement on a pure passive probe.  Given
that root, both outcome-dependent children are optimized exactly by the
one-batch Schur-multiplier SDP.  Root probe and measurement vectors are then
updated by alternating top-eigenvector steps.  This is a falsification
search over a strict subset of all adaptive protocols, not a theorem proof.
"""

from __future__ import annotations

import argparse
import itertools
import json
import platform

import cvxpy as cp
import numpy as np


Q = 2
N = Q * Q
MODES = 4 * N
SEED = 2026071410


def sylvester_sign(n):
    h = np.array([[1]], dtype=int)
    while h.shape[0] < n:
        h = np.block([[h, h], [h, -h]])
    return h


def exact_hard_pair():
    k = sylvester_sign(Q)
    pairs = []
    for permutation in itertools.permutations(range(Q)):
        for signs in itertools.product((-1, 1), repeat=Q):
            p = np.zeros((Q, Q), dtype=int)
            for column, row in enumerate(permutation):
                p[row, column] = signs[column]
            x = k @ p
            y = p @ k
            if not (np.all(np.abs(x) == 1) and np.all(np.abs(y) == 1)):
                raise AssertionError("signed-permutation pair is not Boolean")
            pairs.append((x.reshape(-1), y.reshape(-1)))

    plus_counts = {}
    minus_counts = {}
    for first, second, third in itertools.product(pairs, repeat=3):
        x1, y1 = first
        x2, y2 = second
        x3, y3 = third
        plus = tuple(np.concatenate((x1, y1 * x2, y2 * x3, y3)).tolist())
        minus = tuple(np.concatenate((-x1, y1 * x2, y2 * x3, y3)).tolist())
        plus_counts[plus] = plus_counts.get(plus, 0) + 1
        minus_counts[minus] = minus_counts.get(minus, 0) + 1

    universe = sorted(set(plus_counts) | set(minus_counts))
    normalization = len(pairs) ** 3
    inputs = np.array(universe, dtype=float)
    signed_prior = np.array(
        [(plus_counts.get(x, 0) - minus_counts.get(x, 0)) / normalization for x in universe],
        dtype=float,
    )
    if not np.isclose(signed_prior.sum(), 0.0, atol=1e-14):
        raise AssertionError("signed prior has nonzero mass")
    return inputs, signed_prior


def parity_supports(dose):
    return [
        support
        for size in range(dose + 1)
        for support in itertools.combinations(range(MODES), size)
    ]


def character_matrix(inputs, dose):
    supports = parity_supports(dose)
    result = np.ones((len(inputs), len(supports)), dtype=float)
    for column, support in enumerate(supports):
        if support:
            result[:, column] = np.prod(inputs[:, support], axis=1)
    return result


class OneBatchSDP:
    def __init__(self, characters):
        self.characters = characters
        dimension = characters.shape[1]
        self.symbol = cp.Parameter((dimension, dimension), symmetric=True)
        self.q = cp.Variable(dimension, nonneg=True)
        self.x = cp.Variable((dimension, dimension), symmetric=True)
        diagonal = cp.diag(self.q)
        self.lower = diagonal - self.x >> 0
        self.upper = diagonal + self.x >> 0
        self.problem = cp.Problem(
            cp.Maximize(cp.sum(cp.multiply(self.symbol, self.x))),
            [cp.sum(self.q) == 1, self.lower, self.upper],
        )

    def solve(self, signed_weights):
        c = self.characters.T @ (signed_weights[:, None] * self.characters)
        c = (c + c.T) / 2
        self.symbol.value = c
        optimum = self.problem.solve(
            solver="CLARABEL",
            tol_gap_abs=5e-8,
            tol_gap_rel=5e-8,
            tol_feas=5e-8,
            warm_start=True,
        )
        if self.problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
            raise RuntimeError(self.problem.status)

        # CLARABEL's dual matrices satisfy A-B=C.  Turn their numerical
        # values into an explicitly feasible dual pair by PSD shifts and an
        # eigen-decomposition of the equality residual.  The resulting t is
        # a safe floating-point upper certificate for this fixed branch.
        a = np.array(self.lower.dual_value, dtype=float)
        b = np.array(self.upper.dual_value, dtype=float)
        a = (a + a.T) / 2
        b = (b + b.T) / 2
        safety = 2e-8
        a += max(0.0, -float(np.linalg.eigvalsh(a)[0])) * np.eye(len(a))
        b += max(0.0, -float(np.linalg.eigvalsh(b)[0])) * np.eye(len(b))
        residual = c - (a - b)
        evals, evecs = np.linalg.eigh((residual + residual.T) / 2)
        positive = np.maximum(evals, 0)
        negative = np.maximum(-evals, 0)
        a += (evecs * positive) @ evecs.T
        b += (evecs * negative) @ evecs.T
        dual_upper = float(np.max(np.diag(a + b))) + safety
        primal_tv = max(0.0, float(optimum) / 2)
        dual_tv = max(primal_tv, dual_upper / 2)
        return (
            primal_tv,
            dual_tv,
            np.array(self.x.value, dtype=float),
            np.array(self.q.value, dtype=float),
        )


def outcome_probability(characters, probe, receiver):
    amplitudes = (characters * probe[None, :]) @ np.conj(receiver)
    result = np.abs(amplitudes) ** 2
    if result.min() < -1e-10 or result.max() > 1 + 1e-8:
        raise AssertionError((result.min(), result.max()))
    return np.clip(result, 0, 1)


def top_eigenvector(weight, characters, fixed_vector):
    phased = characters * fixed_vector[None, :]
    matrix = np.einsum("x,xi,xj->ij", weight, phased, phased.conj(), optimize=True)
    matrix = (matrix + matrix.conj().T) / 2
    _, vectors = np.linalg.eigh(matrix)
    return vectors[:, -1]


def one_start(rng, root_characters, child_characters, signed_prior, rounds):
    dimension = root_characters.shape[1]
    probe = rng.normal(size=dimension) + 1j * rng.normal(size=dimension)
    receiver = rng.normal(size=dimension) + 1j * rng.normal(size=dimension)
    probe /= np.linalg.norm(probe)
    receiver /= np.linalg.norm(receiver)
    child_sdp = OneBatchSDP(child_characters)
    best = None

    for _ in range(rounds):
        p_one = outcome_probability(root_characters, probe, receiver)
        lower_one, upper_one, x_one, q_one = child_sdp.solve(signed_prior * p_one)
        lower_zero, upper_zero, x_zero, q_zero = child_sdp.solve(signed_prior * (1 - p_one))
        candidate = {
            "tv_lower": lower_one + lower_zero,
            "tv_fixed_root_upper": upper_one + upper_zero,
            "probe": probe.copy(),
            "receiver": receiver.copy(),
            "q_one": q_one,
            "q_zero": q_zero,
        }
        if best is None or candidate["tv_lower"] > best["tv_lower"]:
            best = candidate

        g_one = np.einsum(
            "xi,ij,xj->x", child_characters, x_one, child_characters, optimize=True
        )
        g_zero = np.einsum(
            "xi,ij,xj->x", child_characters, x_zero, child_characters, optimize=True
        )
        root_weight = 0.5 * signed_prior * (g_one - g_zero)
        for _ in range(4):
            receiver = top_eigenvector(root_weight, root_characters, probe)
            probe = top_eigenvector(root_weight, root_characters, receiver)

    # Score the final root rather than only the pre-update iterates.
    p_one = outcome_probability(root_characters, probe, receiver)
    lower_one, upper_one, _, q_one = child_sdp.solve(signed_prior * p_one)
    lower_zero, upper_zero, _, q_zero = child_sdp.solve(signed_prior * (1 - p_one))
    candidate = {
        "tv_lower": lower_one + lower_zero,
        "tv_fixed_root_upper": upper_one + upper_zero,
        "probe": probe.copy(),
        "receiver": receiver.copy(),
        "q_one": q_one,
        "q_zero": q_zero,
    }
    if candidate["tv_lower"] > best["tv_lower"]:
        best = candidate
    return best


def summarize_amplitudes(vector, supports, count=8):
    order = np.argsort(np.abs(vector))[::-1][:count]
    return [
        {
            "support": list(supports[index]),
            "magnitude": float(abs(vector[index])),
            "phase": float(np.angle(vector[index])),
        }
        for index in order
        if abs(vector[index]) > 1e-7
    ]


def summarize_probabilities(vector, supports, count=8):
    order = np.argsort(vector)[::-1][:count]
    return [
        {"support": list(supports[index]), "weight": float(vector[index])}
        for index in order
        if vector[index] > 1e-7
    ]


def search_split(inputs, signed_prior, root_dose, child_dose, starts, rounds, seed):
    root_supports = parity_supports(root_dose)
    child_supports = parity_supports(child_dose)
    root_characters = character_matrix(inputs, root_dose)
    child_characters = character_matrix(inputs, child_dose)
    rng = np.random.default_rng(seed)
    records = []
    for _ in range(starts):
        records.append(
            one_start(rng, root_characters, child_characters, signed_prior, rounds)
        )
    best = max(records, key=lambda record: record["tv_lower"])
    return {
        "split": f"{root_dose}+{child_dose}",
        "root_dimension": root_characters.shape[1],
        "child_dimension": child_characters.shape[1],
        "starts": starts,
        "rounds": rounds,
        "best_tv_lower": best["tv_lower"],
        "best_fixed_root_dual_upper": best["tv_fixed_root_upper"],
        "equal_prior_error_upper": (1 - best["tv_lower"]) / 2,
        "root_probe_largest_amplitudes": summarize_amplitudes(best["probe"], root_supports),
        "root_effect_largest_amplitudes": summarize_amplitudes(best["receiver"], root_supports),
        "child_q_outcome_one": summarize_probabilities(best["q_one"], child_supports),
        "child_q_outcome_zero": summarize_probabilities(best["q_zero"], child_supports),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--starts", type=int, default=3)
    parser.add_argument("--rounds", type=int, default=4)
    args = parser.parse_args()
    inputs, signed_prior = exact_hard_pair()

    # Sanity certificate: the q=2 crossing-pair construction has value 1/q.
    dose_two = OneBatchSDP(character_matrix(inputs, 2))
    dose_two_lower, dose_two_upper, _, _ = dose_two.solve(signed_prior)
    if not (abs(dose_two_lower - 1 / Q) < 2e-6 and dose_two_upper < 1 / Q + 2e-6):
        raise AssertionError((dose_two_lower, dose_two_upper))

    results = [
        search_split(inputs, signed_prior, 1, 1, args.starts, args.rounds, SEED + 2),
    ]
    report = {
        "scope": "exact q=2 (N=4) plant; binary pure root and rank-one effect; exact SDP children",
        "seed": SEED,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "cvxpy": cp.__version__,
        "solvers": cp.installed_solvers(),
        "dose_two_sanity": {
            "primal_lower": dose_two_lower,
            "repaired_dual_upper": dose_two_upper,
        },
        "searches": results,
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
