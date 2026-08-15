#!/usr/bin/env python3
"""Random stress test of the collision-aware bilateral frame lemma."""

from __future__ import annotations

import numpy as np


SEED = 2026071404
TOL = 3e-9


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
    return (eigenvectors * np.sqrt(np.maximum(eigenvalues, 0.0))) @ eigenvectors.conj().T


def fiber_size(mapping: np.ndarray, omega: int) -> int:
    largest = 0
    for outer in range(mapping.shape[0]):
        counts = np.bincount(mapping[outer], minlength=omega)
        largest = max(largest, int(counts.max()))
    return largest


def stress(trials: int = 400) -> float:
    rng = np.random.default_rng(SEED)
    worst = 0.0

    for _ in range(trials):
        omega = int(rng.integers(1, 7))
        r_count = int(rng.integers(1, 4))
        s_count = int(rng.integers(1, 4))
        xdim = int(rng.integers(1, 6))
        zdim = int(rng.integers(1, 6))
        tau = rng.integers(0, omega, size=(r_count, xdim))
        upsilon = rng.integers(0, omega, size=(s_count, zdim))

        qdiag = rng.random(omega)
        qdiag /= qdiag.sum()
        frame = np.sqrt(qdiag)[:, None] * random_unitary(rng, omega)
        ps = [random_psd(rng, xdim) for _ in range(s_count)]
        qs = [random_psd(rng, zdim) for _ in range(r_count)]

        a_diag = 0.0
        for r in range(r_count):
            for s in range(s_count):
                for x in range(xdim):
                    a_diag += ps[s][x, x].real * qdiag[tau[r, x]]
        b_diag = 0.0
        for s in range(s_count):
            for r in range(r_count):
                for z in range(zdim):
                    b_diag += qs[r][z, z].real * qdiag[upsilon[s, z]]

        kappa_tau = fiber_size(tau, omega)
        kappa_upsilon = fiber_size(upsilon, omega)
        bound = np.sqrt(kappa_tau * kappa_upsilon * a_diag * b_diag)
        lhs = 0.0

        for y in range(omega):
            for r in range(r_count):
                u = frame[tau[r], y]
                for s in range(s_count):
                    v = np.conj(frame[upsilon[s], y])
                    atom = np.kron(u, v)
                    majorant = np.kron(ps[s], qs[r])
                    direction = rng.normal(size=majorant.shape[0])
                    direction = direction + 1j * rng.normal(size=direction.shape)
                    direction /= np.linalg.norm(direction)
                    w = psd_sqrt(majorant) @ direction
                    lhs += abs(np.vdot(atom, w))

        ratio = lhs / bound if bound else 0.0
        worst = max(worst, ratio)
        if ratio > 1 + 10 * TOL:
            raise AssertionError(("collision-aware lemma", ratio, kappa_tau, kappa_upsilon))

    return worst


def main() -> None:
    worst = stress()
    print(f"collision-aware bilateral stress passed: worst_ratio={worst:.12g}")


if __name__ == "__main__":
    main()
