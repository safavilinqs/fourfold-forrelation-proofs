#!/usr/bin/env python3
"""Random and saturating checks of the universal weighted link Gram bound."""

from __future__ import annotations

import numpy as np


SEED = 2026071421


def main() -> None:
    rng = np.random.default_rng(SEED)
    worst_ratio = 0.0
    checked = 0
    for latent in (1, 2, 5, 11, 32):
        for left in (1, 3, 9, 25):
            for right in (1, 4, 10, 23):
                for _ in range(8):
                    probability = rng.dirichlet(np.ones(latent))
                    left_features = np.exp(
                        2j * np.pi * rng.random((latent, left))
                    )
                    right_features = np.exp(
                        2j * np.pi * rng.random((latent, right))
                    )
                    moment = (
                        left_features.conj().T
                        @ (probability[:, None] * right_features)
                    )
                    left_mass = float(rng.random())
                    right_mass = float(rng.random())
                    p = left_mass * rng.dirichlet(np.ones(left))
                    q = right_mass * rng.dirichlet(np.ones(right))
                    weighted = (
                        np.sqrt(p)[:, None] * moment * np.sqrt(q)[None, :]
                    )
                    nuclear = float(
                        np.linalg.svd(weighted, compute_uv=False).sum()
                    )
                    bound = np.sqrt(left_mass * right_mass)
                    if nuclear > bound * (1 + 4e-12):
                        raise AssertionError(
                            (
                                "weighted link Gram contraction",
                                latent,
                                left,
                                right,
                                nuclear,
                                bound,
                            )
                        )
                    if bound:
                        worst_ratio = max(worst_ratio, nuclear / bound)
                    checked += 1

    # One latent atom makes every feature column collinear and saturates
    # Schatten Hölder for point-mass diagonal weights.
    moment = np.ones((3, 4))
    p = np.array([1.0, 0.0, 0.0])
    q = np.array([1.0, 0.0, 0.0, 0.0])
    saturation = float(
        np.linalg.svd(np.sqrt(p)[:, None] * moment * np.sqrt(q)[None, :],
                      compute_uv=False).sum()
    )
    if saturation != 1:
        raise AssertionError(("saturation", saturation))

    print(
        "weighted link Gram contraction passed: "
        f"instances={checked}, worst_random_ratio={worst_ratio:.12g}, "
        f"saturation={saturation:.12g}"
    )


if __name__ == "__main__":
    main()
