#!/usr/bin/env python3
"""Exact regression for the Hilbert-versus-projective ancillary norm gap."""

import numpy as np


def sylvester(n: int) -> np.ndarray:
    h = np.array([[1.0]])
    while h.shape[0] < n:
        h = np.block([[h, h], [h, -h]])
    return h / np.sqrt(n)


def main() -> None:
    results = []
    for n in (2, 4, 8, 16):
        h = sylvester(n)
        injective = np.linalg.norm(h, ord=2)
        auxiliary = h.reshape(-1) / np.sqrt(n)
        hilbert = np.linalg.norm(auxiliary)
        contraction = abs(np.vdot(h.reshape(-1), auxiliary))
        projective = np.linalg.norm(auxiliary.reshape(n, n), ord="nuc")
        if not np.isclose(injective, 1.0, atol=2e-12):
            raise AssertionError(("edge injective norm", n, injective))
        if not np.isclose(hilbert, 1.0, atol=2e-12):
            raise AssertionError(("auxiliary Hilbert norm", n, hilbert))
        if not np.isclose(projective, np.sqrt(n), atol=2e-12):
            raise AssertionError(("auxiliary projective norm", n, projective))
        if not np.isclose(contraction, np.sqrt(n), atol=2e-12):
            raise AssertionError(("dimension amplification", n, contraction))
        results.append(f"N={n}:gap={contraction:.12g}")
    print("singleton ancillary norm gap confirmed: " + ", ".join(results))


if __name__ == "__main__":
    main()
