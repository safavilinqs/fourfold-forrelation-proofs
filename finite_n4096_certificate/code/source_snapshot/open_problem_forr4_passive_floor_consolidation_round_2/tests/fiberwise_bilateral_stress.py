#!/usr/bin/env python3
"""Adversarial checks for the corrected fiberwise bilateral lemma.

The script does two different things:

1. records an exact counterexample to the unqualified arbitrary-map lemma;
2. stress-tests the injective scalar and Hilbert-valued forms with general
   entangled selected tensors satisfying the required operator majorant.

Passing is regression evidence, not a proof of the reverse-tree induction.
"""

from __future__ import annotations

import numpy as np


SEED = 2026071402
TOL = 2e-9


def random_unitary(rng: np.random.Generator, n: int) -> np.ndarray:
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, r = np.linalg.qr(z)
    diagonal = np.diag(r)
    phase = np.where(np.abs(diagonal) > 0, diagonal / np.abs(diagonal), 1.0)
    return q * phase.conj()


def random_psd(rng: np.random.Generator, n: int) -> np.ndarray:
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    value = z @ z.conj().T
    return value / max(np.linalg.norm(value, ord=2), 1e-15)


def psd_sqrt(value: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(value)
    eigenvalues = np.maximum(eigenvalues, 0.0)
    return (eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.conj().T


def collision_counterexample() -> float:
    # U=(1,1), V=(1), P=all-ones, Q=(1), W=(1,1).
    u = np.ones(2, dtype=complex)
    v = np.ones(1, dtype=complex)
    p = np.ones((2, 2), dtype=complex)
    q = np.ones((1, 1), dtype=complex)
    w = np.ones(2, dtype=complex)
    majorant_residual = np.linalg.eigvalsh(np.kron(p, q) - np.outer(w, w.conj())).min()
    if majorant_residual < -TOL:
        raise AssertionError(("collision majorant", majorant_residual))
    lhs = abs(np.vdot(np.kron(u, v), w))
    diagonal_bound = np.sqrt(np.trace(p).real * np.trace(q).real)
    ratio = lhs / diagonal_bound
    if not np.isclose(ratio, np.sqrt(2.0), atol=1e-12):
        raise AssertionError(("unexpected collision ratio", ratio))
    exact_collision_bound = np.sqrt(p.sum().real * q.sum().real)
    if not np.isclose(lhs, exact_collision_bound, atol=1e-12):
        raise AssertionError(("collision-aware bound", lhs, exact_collision_bound))
    return ratio


def injective_maps(rng: np.random.Generator, outer: int, inner: int, omega: int) -> np.ndarray:
    if inner > omega:
        raise ValueError("injective slice is impossible")
    return np.stack([rng.choice(omega, size=inner, replace=False) for _ in range(outer)])


def lambda_bound(tau: np.ndarray, ps: list[np.ndarray], omega: int) -> float:
    values = np.zeros(omega)
    diagonal_sum = sum(np.diag(p).real for p in ps)
    for r in range(tau.shape[0]):
        for x in range(tau.shape[1]):
            values[tau[r, x]] += diagonal_sum[x]
    return float(values.max())


def mu_bound(upsilon: np.ndarray, qs: list[np.ndarray], omega: int) -> float:
    values = np.zeros(omega)
    diagonal_sum = sum(np.diag(q).real for q in qs)
    for s in range(upsilon.shape[0]):
        for z in range(upsilon.shape[1]):
            values[upsilon[s, z]] += diagonal_sum[z]
    return float(values.max())


def stress_injective_forms(trials: int = 300) -> tuple[float, float, float]:
    rng = np.random.default_rng(SEED)
    worst_scalar = 0.0
    worst_vector = 0.0
    worst_majorant_slack = 0.0

    for _ in range(trials):
        omega = int(rng.integers(4, 8))
        r_count = int(rng.integers(1, 4))
        s_count = int(rng.integers(1, 4))
        x_count = int(rng.integers(1, min(4, omega) + 1))
        z_count = int(rng.integers(1, min(4, omega) + 1))
        auxiliary = int(rng.integers(1, 5))

        tau = injective_maps(rng, r_count, x_count, omega)
        upsilon = injective_maps(rng, s_count, z_count, omega)

        qdiag = rng.random(omega)
        qdiag /= qdiag.sum()
        # Columns are frame vectors, so frame @ frame^* = D_q.
        frame = np.sqrt(qdiag)[:, None] * random_unitary(rng, omega)
        if not np.allclose(frame @ frame.conj().T, np.diag(qdiag), atol=2e-12):
            raise AssertionError("frame completeness")

        ps = [random_psd(rng, x_count) for _ in range(s_count)]
        qs = [random_psd(rng, z_count) for _ in range(r_count)]
        bound = np.sqrt(lambda_bound(tau, ps, omega) * mu_bound(upsilon, qs, omega))
        if bound < 1e-15:
            continue

        scalar_lhs = 0.0
        vector_sum = np.zeros(auxiliary, dtype=complex)

        for y in range(omega):
            for r in range(r_count):
                u = frame[tau[r], y]
                for s in range(s_count):
                    v = np.conj(frame[upsilon[s], y])
                    atom = np.kron(u, v)
                    majorant = np.kron(ps[s], qs[r])
                    root = psd_sqrt(majorant)
                    dimension = majorant.shape[0]

                    # A generic entangled scalar selected vector.
                    direction = rng.normal(size=dimension) + 1j * rng.normal(size=dimension)
                    direction /= np.linalg.norm(direction)
                    w = root @ direction
                    residual = majorant - np.outer(w, w.conj())
                    slack = max(0.0, -float(np.linalg.eigvalsh(residual).min()))
                    worst_majorant_slack = max(worst_majorant_slack, slack)
                    if slack > TOL:
                        raise AssertionError(("scalar operator majorant", slack))
                    scalar_lhs += abs(np.vdot(atom, w))

                    # A Hilbert-valued selected tensor W with W W^* <= majorant.
                    zmap = rng.normal(size=(dimension, auxiliary)) + 1j * rng.normal(
                        size=(dimension, auxiliary)
                    )
                    operator_norm = np.linalg.norm(zmap, ord=2)
                    if operator_norm > 0:
                        zmap /= operator_norm
                    wmap = root @ zmap
                    residual = majorant - wmap @ wmap.conj().T
                    slack = max(0.0, -float(np.linalg.eigvalsh(residual).min()))
                    worst_majorant_slack = max(worst_majorant_slack, slack)
                    if slack > 5 * TOL:
                        raise AssertionError(("vector operator majorant", slack))
                    phase = np.exp(2j * np.pi * rng.random())
                    vector_sum += phase * (wmap.conj().T @ atom)

        scalar_ratio = scalar_lhs / bound
        vector_ratio = np.linalg.norm(vector_sum) / bound
        worst_scalar = max(worst_scalar, scalar_ratio)
        worst_vector = max(worst_vector, vector_ratio)
        if scalar_ratio > 1 + 10 * TOL:
            raise AssertionError(("injective scalar lemma", scalar_ratio))
        if vector_ratio > 1 + 10 * TOL:
            raise AssertionError(("injective Hilbert-valued lemma", vector_ratio))

    return worst_scalar, worst_vector, worst_majorant_slack


def main() -> None:
    collision_ratio = collision_counterexample()
    scalar, vector, slack = stress_injective_forms()
    print(
        "fiberwise bilateral stress completed: "
        f"expected_collision_ratio={collision_ratio:.12g}, "
        f"injective_scalar_worst={scalar:.12g}, "
        f"injective_vector_worst={vector:.12g}, "
        f"max_numerical_majorant_slack={slack:.3g}"
    )


if __name__ == "__main__":
    main()
