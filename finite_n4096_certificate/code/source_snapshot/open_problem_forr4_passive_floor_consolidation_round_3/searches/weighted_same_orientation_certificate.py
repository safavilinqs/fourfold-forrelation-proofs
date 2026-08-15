#!/usr/bin/env python3
"""Translation-twirled arbitrary-diagonal certificate for one endpoint slice.

The occurrence orientation has rows (i,b,d) and columns (E,F,c).  Translation
twirling and joint concavity reduce arbitrary diagonal row/column laws to a
joint law Q[x,y] on the two nonzero endpoint-pair XORs.  The complete
spectrum and supporting gradient are computed from exact integer Walsh data.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from alternating_double_endpoint_spectrum import (
    checked_fwht,
    orbit_correlation,
)


@dataclass(frozen=True)
class SameOrientationCertificate:
    order: int
    coefficient: float
    supporting_upper: float
    high_orbit_mass: float
    active_gradient_spread: float
    inactive_gradient_gap: float


def response_data(order: int) -> tuple[np.ndarray, np.ndarray, float]:
    dimension = order * order
    correlation = orbit_correlation(order)
    endpoint_symbol = checked_fwht(correlation, axis=0).astype(float)
    scale = 4 / (dimension**8 * (order - 1) ** 4)
    return correlation, endpoint_symbol, scale


def vertical_orbit_spectra(
    order: int,
    endpoint_symbol: np.ndarray,
    scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    dimension = order * order
    vertical = np.arange(order, dimension, order, dtype=np.int64)
    dot = np.empty((len(vertical), len(vertical)), dtype=np.int8)
    for row, left in enumerate(vertical):
        for column, right in enumerate(vertical):
            dot[row, column] = (
                int(left & right).bit_count() % 2
            )
    high_mask = dot == 0
    low_mask = ~high_mask
    high_law = high_mask / int(high_mask.sum())
    low_law = low_mask / int(low_mask.sum())
    restricted_symbol = endpoint_symbol[:, vertical]
    high_spectrum = np.empty((dimension, dimension))
    low_spectrum = np.empty((dimension, dimension))

    for alpha in range(dimension):
        shifted = endpoint_symbol[
            np.bitwise_xor(alpha, vertical)[:, None],
            vertical[None, :],
        ]
        high_spectrum[alpha] = scale * (
            restricted_symbol @ (high_law * shifted).sum(axis=1)
        )
        low_spectrum[alpha] = scale * (
            restricted_symbol @ (low_law * shifted).sum(axis=1)
        )
    return (
        high_spectrum,
        low_spectrum,
        high_mask,
        low_mask,
    )


def choose_orbit_mass(
    order: int,
    high_spectrum: np.ndarray,
    low_spectrum: np.ndarray,
) -> float:
    dimension = order * order
    high = high_spectrum.reshape(-1)
    low = low_spectrum.reshape(-1)

    def derivative(mass: float) -> float:
        spectrum = mass * high + (1 - mass) * low
        return float(
            dimension
            * np.sum((high - low) / (2 * np.sqrt(spectrum)))
        )

    tolerance = 1e-12
    lower = tolerance
    upper = 1 - tolerance
    if derivative(upper) >= 0:
        return 1.0
    return float(brentq(derivative, lower, upper, xtol=1e-14))


def full_difference_gradient(
    order: int,
    endpoint_symbol: np.ndarray,
    scale: float,
    spectrum: np.ndarray,
) -> np.ndarray:
    dimension = order * order
    reciprocal_root = 1 / np.sqrt(spectrum)
    gradient = np.empty((dimension, dimension))
    factor = dimension * scale / 2
    indices = np.arange(dimension)
    for left_difference in range(dimension):
        contracted = reciprocal_root @ endpoint_symbol[:, left_difference]
        gradient[left_difference] = factor * (
            endpoint_symbol[indices ^ left_difference].T @ contracted
        )
    return gradient


def certificate(order: int) -> SameOrientationCertificate:
    dimension = order * order
    _, endpoint_symbol, scale = response_data(order)
    high, low, high_mask, low_mask = vertical_orbit_spectra(
        order, endpoint_symbol, scale
    )
    mass = choose_orbit_mass(order, high, low)
    spectrum = mass * high + (1 - mass) * low
    if np.any(spectrum <= 0):
        raise AssertionError(("positive weighted spectrum", order))
    coefficient = float(dimension * np.sqrt(spectrum).sum())
    gradient = full_difference_gradient(
        order, endpoint_symbol, scale, spectrum
    )

    vertical = np.arange(order, dimension, order, dtype=np.int64)
    active_pairs: list[tuple[int, int]] = []
    if mass > 1e-10:
        active_pairs.extend(
            (int(vertical[row]), int(vertical[column]))
            for row, column in zip(*np.nonzero(high_mask), strict=True)
        )
    if 1 - mass > 1e-10:
        active_pairs.extend(
            (int(vertical[row]), int(vertical[column]))
            for row, column in zip(*np.nonzero(low_mask), strict=True)
        )
    active_values = np.array(
        [gradient[left, right] for left, right in active_pairs]
    )
    level = coefficient / 2
    active_spread = float(np.ptp(active_values))

    active_set = set(active_pairs)
    inactive_maximum = max(
        gradient[left, right]
        for left in range(1, dimension)
        for right in range(1, dimension)
        if (left, right) not in active_set
    )
    inactive_gap = float(inactive_maximum - level)
    maximum_gradient = float(
        np.max(gradient[1:, 1:])
    )

    if active_spread > 2e-12:
        raise AssertionError(("active KKT spread", order, active_spread))
    if inactive_gap > 2e-12:
        raise AssertionError(("inactive KKT violation", order, inactive_gap))

    # For a concave degree-one-half objective in Q, the tangent upper is
    # Phi/2 + max grad.  Retain a conservative numerical allowance.
    supporting_upper = coefficient / 2 + maximum_gradient + 1e-12
    return SameOrientationCertificate(
        order=order,
        coefficient=coefficient,
        supporting_upper=supporting_upper,
        high_orbit_mass=mass,
        active_gradient_spread=active_spread,
        inactive_gradient_gap=inactive_gap,
    )


def main() -> None:
    for order in (4, 8, 16, 32):
        result = certificate(order)
        print(
            f"q={order},N={order * order},"
            f"high_orbit_mass={result.high_orbit_mass:.15g},"
            f"coefficient={result.coefficient:.15g},"
            f"supporting_upper={result.supporting_upper:.15g},"
            f"active_spread={result.active_gradient_spread:.3g},"
            f"inactive_gap={result.inactive_gradient_gap:.12g}"
        )


if __name__ == "__main__":
    main()
