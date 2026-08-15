#!/usr/bin/env python3
"""Two-level adaptive stress test for the reverse preimage invariant.

This is an abstract operator test of the exact step used by the proof.  A
root outcome selects a different complete child frame.  The child removes a
bilateral marked frame while retaining an auxiliary output in the range of
the root graph map.  The test checks both the returned preimage norm and the
final root contraction.

It is falsification evidence only: the physical skeleton still has to be
shown to instantiate these operator hypotheses.
"""

from __future__ import annotations

import numpy as np


SEED = 2026071403
TOL = 5e-9


def random_unitary(rng: np.random.Generator, n: int) -> np.ndarray:
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, r = np.linalg.qr(z)
    diagonal = np.diag(r)
    phase = np.where(np.abs(diagonal) > 0, diagonal / np.abs(diagonal), 1.0)
    return q * phase.conj()


def complete_frame(rng: np.random.Generator, omega: int) -> np.ndarray:
    qdiag = rng.random(omega)
    qdiag /= qdiag.sum()
    frame = np.sqrt(qdiag)[:, None] * random_unitary(rng, omega)
    if not np.allclose(frame @ frame.conj().T, np.diag(qdiag), atol=2e-12):
        raise AssertionError("frame completeness")
    return frame


def random_psd(rng: np.random.Generator, n: int) -> np.ndarray:
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    value = z @ z.conj().T
    return value / max(np.linalg.norm(value, ord=2), 1e-15)


def psd_sqrt(value: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(value)
    return (eigenvectors * np.sqrt(np.maximum(eigenvalues, 0.0))) @ eigenvectors.conj().T


def injective_maps(rng: np.random.Generator, outer: int, inner: int, omega: int) -> np.ndarray:
    return np.stack([rng.choice(omega, size=inner, replace=False) for _ in range(outer)])


def local_bounds(
    tau: np.ndarray,
    upsilon: np.ndarray,
    ps: list[np.ndarray],
    qs: list[np.ndarray],
    omega: int,
) -> float:
    lambda_values = np.zeros(omega)
    pdiag = sum(np.diag(p).real for p in ps)
    for r in range(tau.shape[0]):
        for x in range(tau.shape[1]):
            lambda_values[tau[r, x]] += pdiag[x]

    mu_values = np.zeros(omega)
    qdiag = sum(np.diag(q).real for q in qs)
    for s in range(upsilon.shape[0]):
        for z in range(upsilon.shape[1]):
            mu_values[upsilon[s, z]] += qdiag[z]
    return float(np.sqrt(lambda_values.max() * mu_values.max()))


def child_elimination(
    rng: np.random.Generator,
    ancestor_map: np.ndarray,
    omega: int,
    xdim: int,
    zdim: int,
    r_count: int,
    s_count: int,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Return physical output, its explicit preimage, bound, and majorant slack."""

    frame = complete_frame(rng, omega)
    tau = injective_maps(rng, r_count, xdim, omega)
    upsilon = injective_maps(rng, s_count, zdim, omega)
    ps = [random_psd(rng, xdim) for _ in range(s_count)]
    qs = [random_psd(rng, zdim) for _ in range(r_count)]
    bound = local_bounds(tau, upsilon, ps, qs, omega)
    preimage = np.zeros(ancestor_map.shape[1], dtype=complex)
    physical = np.zeros(ancestor_map.shape[0], dtype=complex)
    worst_slack = 0.0

    for y in range(omega):
        for r in range(r_count):
            u = frame[tau[r], y]
            for s in range(s_count):
                v = np.conj(frame[upsilon[s], y])
                atom = np.kron(u, v)
                majorant = np.kron(ps[s], qs[r])
                root = psd_sqrt(majorant)
                zmap = rng.normal(size=(majorant.shape[0], ancestor_map.shape[1]))
                zmap = zmap + 1j * rng.normal(size=zmap.shape)
                norm = np.linalg.norm(zmap, ord=2)
                if norm:
                    zmap /= norm

                # W is current-selected coordinates by physical ancestor output.
                wmap = root @ zmap @ ancestor_map.conj().T
                residual = majorant - wmap @ wmap.conj().T
                slack = max(0.0, -float(np.linalg.eigvalsh(residual).min()))
                worst_slack = max(worst_slack, slack)
                if slack > TOL:
                    raise AssertionError(("child majorant", slack))

                phase = np.exp(2j * np.pi * rng.random())
                increment = phase * (zmap.conj().T @ root @ atom)
                preimage += increment
                physical += phase * (wmap.conj().T @ atom)

    if not np.allclose(physical, ancestor_map @ preimage, atol=2e-10):
        raise AssertionError("child output left ancestor range")
    if np.linalg.norm(preimage) > bound * (1 + 10 * TOL):
        raise AssertionError(("returned preimage norm", np.linalg.norm(preimage), bound))
    return physical, preimage, bound, worst_slack


def stress_two_level(trials: int = 250) -> tuple[float, float, float]:
    rng = np.random.default_rng(SEED)
    worst_preimage = 0.0
    worst_tree = 0.0
    worst_slack = 0.0

    for _ in range(trials):
        selected_dim = int(rng.integers(2, 5))
        preimage_dim = int(rng.integers(1, selected_dim + 1))
        ancestor_map = rng.normal(size=(selected_dim, preimage_dim))
        ancestor_map = ancestor_map + 1j * rng.normal(size=ancestor_map.shape)
        ancestor_map /= max(np.linalg.norm(ancestor_map, ord=2), 1e-15)
        ancestor_covariance = ancestor_map @ ancestor_map.conj().T

        root_frame = complete_frame(rng, selected_dim + 1)
        child_outputs = []
        child_bounds = []

        # The child instrument is chosen adaptively from the root outcome.
        for root_outcome in range(selected_dim + 1):
            omega = int(rng.integers(4, 8))
            xdim = int(rng.integers(1, min(3, omega) + 1))
            zdim = int(rng.integers(1, min(3, omega) + 1))
            r_count = int(rng.integers(1, 3))
            s_count = int(rng.integers(1, 3))
            output, preimage, bound, slack = child_elimination(
                rng, ancestor_map, omega, xdim, zdim, r_count, s_count
            )
            child_outputs.append(output)
            child_bounds.append(bound)
            worst_slack = max(worst_slack, slack)
            if bound:
                worst_preimage = max(worst_preimage, np.linalg.norm(preimage) / bound)

        child_bound = max(child_bounds)
        root_lambda = float(np.max(np.diag(ancestor_covariance).real))
        total = 0.0j
        absolute_total = 0.0
        for outcome, output in enumerate(child_outputs):
            base = root_frame[0, outcome]
            selected = root_frame[1:, outcome]
            term = np.conj(base) * np.vdot(selected, output)
            total += term
            absolute_total += abs(term)

        root_bound = child_bound * np.sqrt(root_lambda)
        if root_bound:
            ratio = absolute_total / root_bound
            worst_tree = max(worst_tree, ratio)
            if ratio > 1 + 20 * TOL:
                raise AssertionError(("two-level reverse contraction", ratio, abs(total)))

    return worst_preimage, worst_tree, worst_slack


def main() -> None:
    preimage, tree, slack = stress_two_level()
    print(
        "reverse-tree operator stress passed: "
        f"preimage_worst={preimage:.12g}, "
        f"two_level_worst={tree:.12g}, "
        f"max_numerical_majorant_slack={slack:.3g}"
    )


if __name__ == "__main__":
    main()
