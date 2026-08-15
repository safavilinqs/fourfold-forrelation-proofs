#!/usr/bin/env python3
"""Stress the operator-valued form of the reverse-frontier invariant."""

from __future__ import annotations

import numpy as np


SEED = 2026071412
TOL = 3e-9


def random_orthogonal(rng: np.random.Generator, n: int) -> np.ndarray:
    q, r = np.linalg.qr(rng.normal(size=(n, n)))
    signs = np.where(np.diag(r) >= 0, 1.0, -1.0)
    return q * signs


def random_psd(rng: np.random.Generator, n: int) -> np.ndarray:
    value = rng.normal(size=(n, n))
    value = value @ value.T
    return value / max(np.linalg.norm(value, ord=2), 1e-15)


def psd_sqrt(value: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(value)
    return (eigenvectors * np.sqrt(np.maximum(eigenvalues, 0.0))) @ eigenvectors.T


def identity_wire_gap() -> list[tuple[int, float, float]]:
    result = []
    for n in (2, 4, 8, 16):
        identity = np.eye(n)
        operator = float(np.linalg.norm(identity, ord=2))
        hilbert_schmidt = float(np.linalg.norm(identity))
        if not np.isclose(operator, 1.0, atol=1e-14):
            raise AssertionError(("identity operator norm", n, operator))
        if not np.isclose(hilbert_schmidt, np.sqrt(n), atol=1e-14):
            raise AssertionError(("identity HS norm", n, hilbert_schmidt))
        result.append((n, operator, hilbert_schmidt))
    return result


def stress_operator_update(trials: int = 350) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    worst_ratio = 0.0
    worst_majorant_slack = 0.0

    for _ in range(trials):
        omega = int(rng.integers(3, 8))
        xdim = int(rng.integers(1, omega + 1))
        zdim = int(rng.integers(1, omega + 1))
        retained = int(rng.integers(1, 5))
        open_input = int(rng.integers(1, 5))

        qdiag = rng.random(omega)
        qdiag /= qdiag.sum()
        frame = np.sqrt(qdiag)[:, None] * random_orthogonal(rng, omega)
        tau = rng.choice(omega, size=xdim, replace=False)
        upsilon = rng.choice(omega, size=zdim, replace=False)
        p = random_psd(rng, xdim)
        q = random_psd(rng, zdim)
        majorant = np.kron(p, q)
        root = psd_sqrt(majorant)

        a_diag = float(np.dot(np.diag(p), qdiag[tau]))
        b_diag = float(np.dot(np.diag(q), qdiag[upsilon]))
        bound = np.sqrt(max(0.0, a_diag * b_diag))
        reverse = np.zeros((retained, open_input))

        for outcome in range(omega):
            atom = np.kron(frame[tau, outcome], frame[upsilon, outcome])
            preimage = rng.normal(
                size=(majorant.shape[0] * retained, open_input)
            )
            norm = np.linalg.norm(preimage, ord=2)
            if norm:
                preimage /= norm
            phase = -1.0 if rng.integers(2) else 1.0

            for column in range(open_input):
                vector = preimage[:, column]
                selected = root @ vector.reshape(majorant.shape[0], retained)
                residual = majorant - selected @ selected.T
                slack = max(0.0, -float(np.linalg.eigvalsh(residual).min()))
                worst_majorant_slack = max(worst_majorant_slack, slack)
                if slack > TOL:
                    raise AssertionError(("fixed-input covariance majorant", slack))
                reverse[:, column] += phase * (selected.T @ atom)

        ratio = float(np.linalg.norm(reverse, ord=2) / bound) if bound else 0.0
        worst_ratio = max(worst_ratio, ratio)
        if ratio > 1 + 20 * TOL:
            raise AssertionError(("operator-valued reverse update", ratio, bound))

    return worst_ratio, worst_majorant_slack


def main() -> None:
    gaps = identity_wire_gap()
    ratio, slack = stress_operator_update()
    formatted = ", ".join(
        f"N={n}:op={operator:.0f},HS={hilbert_schmidt:.12g}"
        for n, operator, hilbert_schmidt in gaps
    )
    print(
        "open-frontier operator invariant passed: "
        f"{formatted}, worst_update_ratio={ratio:.12g}, "
        f"max_numerical_majorant_slack={slack:.3g}"
    )


if __name__ == "__main__":
    main()
