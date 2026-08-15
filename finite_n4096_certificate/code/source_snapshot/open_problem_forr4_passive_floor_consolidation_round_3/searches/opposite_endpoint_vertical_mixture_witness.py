#!/usr/bin/env python3
"""Explicit mixed-orbit obstruction for the opposite endpoint chain.

The row law is uniform over all nonzero vertical differences for both
selected endpoint pairs.  The column law is uniform over translation orbits
of triples contained in one hidden column.  A frequency-type reduction
evaluates the exact critical-split coefficient through q=32 using only
fourteen small nuclear norms.

At q=32 this physical law exceeded the coefficient compatible with a
historical independent profile/split Perron ledger.  The later chained
accepted-sector repair gives ample slack even after charging these cuts.
The physical witness remains a permanent falsification of a universal
``1/q`` upper bound and of the old scalar budget; the current output reports
its value under the repaired ledger.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

import numpy as np
from scipy.optimize import minimize_scalar

from attenuation_promise_concentration import (
    euclidean_promise_concentration,
    promise_concentration,
)
from occupation_compatible_sector_optimization import certificate
from opposite_endpoint_orbit_scan import (
    cubic_response,
    quintic_response,
    support_xor,
    triple_orbit_representatives,
    walsh_transform,
)


@dataclass(frozen=True)
class VerticalMixtureWitness:
    order: int
    coefficient: float
    frequency_classes: int


@dataclass(frozen=True)
class ForcedLedgerObstruction:
    coefficient: float
    optimal_beta: float
    optimized_total: float
    threshold_overshoot: float


@dataclass(frozen=True)
class ForcedLedgerRepair:
    coefficient: float
    optimal_beta: float
    optimized_total: float
    threshold_slack: float


def character(left: int, right: int) -> int:
    return -1 if int(left & right).bit_count() % 2 else 1


def vertical_mixture_witness(order: int) -> VerticalMixtureWitness:
    if order < 4 or order & (order - 1):
        raise ValueError(("power-of-two order at least four required", order))
    q = order
    dimension = q * q
    vertical = tuple(row * q for row in range(1, q))
    triples = tuple(
        tuple(row * q for row in triple)
        for triple in triple_orbit_representatives(q)
    )
    cubic = np.asarray(
        [
            walsh_transform(cubic_response(q, difference, False))
            for difference in vertical
        ]
    )
    quintic = np.empty((q - 1, len(triples), dimension))
    frequencies = range(dimension)
    for difference_index, difference in enumerate(vertical):
        for triple_index, triple in enumerate(triples):
            triple_xor = support_xor(triple)
            quintic[difference_index, triple_index] = walsh_transform(
                quintic_response(q, difference, triple, False)
            ) * np.asarray(
                [character(triple_xor, frequency) for frequency in frequencies]
            )

    # Under the row-label linear group, a nonzero alpha-row frequency and
    # nonzero gamma-row frequency are classified by their bilinear pairing.
    # The representatives 1,2 have pairing zero and 1,1 have pairing one.
    row_frequency_classes = (
        (0, 0, 1),
        (0, 1, q - 1),
        (1, 0, q - 1),
        (1, 2, (q - 1) * (q // 2 - 1)),
        (1, 1, (q - 1) * (q // 2)),
    )
    nuclear_sum = 0.0
    class_count = 0
    for alpha_row, gamma_row, row_multiplicity in row_frequency_classes:
        alpha_columns = (
            ((0, 1), (1, q - 1))
            if alpha_row == 0
            else ((0, q),)
        )
        for alpha_column, alpha_multiplicity in alpha_columns:
            for gamma_column, gamma_multiplicity in ((0, 1), (1, q - 1)):
                alpha = alpha_row * q + alpha_column
                gamma = gamma_row * q + gamma_column
                shifted_frequencies = np.asarray(
                    [gamma ^ difference for difference in vertical]
                )
                endpoint = quintic[:, :, shifted_frequencies].transpose(
                    2, 0, 1
                )
                matrix = (
                    cubic[:, alpha, None, None] * endpoint
                ).reshape((q - 1) ** 2, len(triples))
                matrix /= (q - 1) * sqrt(len(triples))
                multiplicity = (
                    row_multiplicity
                    * alpha_multiplicity
                    * gamma_multiplicity
                )
                nuclear_sum += multiplicity * float(
                    np.linalg.svd(matrix, compute_uv=False).sum()
                )
                class_count += 1
    return VerticalMixtureWitness(
        order=q,
        coefficient=nuclear_sum / dimension**3,
        frequency_classes=class_count,
    )


def forced_split_coefficients(
    value: float,
) -> dict[tuple[tuple[int, ...], tuple[int, ...]], float]:
    profile = (3, 1, 1, 5)
    reverse = tuple(reversed(profile))
    split = (2, 0, 1, 2)
    complement = tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )
    return {
        (profile, split): value,
        (profile, complement): value,
        (reverse, tuple(reversed(split))): value,
        (reverse, tuple(reversed(complement))): value,
    }


def forced_ledger_obstruction(
    value: float,
) -> ForcedLedgerObstruction:
    """Evaluate the old forced-cut diagnostic using the current ledger.

    The name is retained for reproducibility.  After the chained accepted-
    sector repair, ``threshold_overshoot`` is negative.
    """
    coefficients = forced_split_coefficients(value)

    def total(beta: float) -> float:
        ledger = certificate(
            beta=beta,
            profile_split_coefficients=coefficients,
        )
        return (
            ledger.supporting_upper
            + promise_concentration(1024, beta).two_hypothesis_loss
        )

    optimum = minimize_scalar(
        total,
        bounds=(0.75, 0.81),
        method="bounded",
        options={"xatol": 1e-12},
    )
    optimized_total = float(optimum.fun)
    return ForcedLedgerObstruction(
        coefficient=value,
        optimal_beta=float(optimum.x),
        optimized_total=optimized_total,
        threshold_overshoot=optimized_total - 1 / 3,
    )


def forced_ledger_euclidean_repair(value: float) -> ForcedLedgerRepair:
    """Reoptimize the forced-cut diagnostic with finite-tilt concentration."""

    coefficients = forced_split_coefficients(value)

    def total(beta: float) -> float:
        ledger = certificate(
            beta=beta,
            profile_split_coefficients=coefficients,
        )
        return (
            ledger.supporting_upper
            + euclidean_promise_concentration(
                1024, beta
            ).two_hypothesis_loss
        )

    optimum = minimize_scalar(
        total,
        bounds=(0.75, 0.81),
        method="bounded",
        options={"xatol": 1e-12},
    )
    optimized_total = float(optimum.fun)
    return ForcedLedgerRepair(
        coefficient=value,
        optimal_beta=float(optimum.x),
        optimized_total=optimized_total,
        threshold_slack=1 / 3 - optimized_total,
    )


def main() -> None:
    witnesses = tuple(
        vertical_mixture_witness(order) for order in (4, 8, 16, 32)
    )
    for witness in witnesses:
        print(
            "opposite endpoint vertical mixture: "
            f"q={witness.order},N={witness.order**2},"
            f"frequency_classes={witness.frequency_classes},"
            f"coefficient={witness.coefficient:.15g}"
        )
    obstruction = forced_ledger_obstruction(witnesses[-1].coefficient)
    print(
        "forced critical-cut diagnostic after accepted-sector repair: "
        f"coefficient={obstruction.coefficient:.15g},"
        f"optimal_beta={obstruction.optimal_beta:.15g},"
        f"optimized_total={obstruction.optimized_total:.15g},"
        f"threshold_slack={-obstruction.threshold_overshoot:.15g}"
    )
    repair = forced_ledger_euclidean_repair(witnesses[-1].coefficient)
    print(
        "forced critical-cut Euclidean promise repair: "
        f"coefficient={repair.coefficient:.15g},"
        f"optimal_beta={repair.optimal_beta:.15g},"
        f"optimized_total={repair.optimized_total:.15g},"
        f"threshold_slack={repair.threshold_slack:.15g}"
    )


if __name__ == "__main__":
    main()
