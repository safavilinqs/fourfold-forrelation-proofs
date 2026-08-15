#!/usr/bin/env python3
"""Quantify the dose-six barrier from summing fixed-split bounds."""

from __future__ import annotations

from itertools import product
from math import comb, sqrt


Q = 32
BETA = 5 / 6


def endpoint_energies() -> list[float]:
    q = Q
    return [
        (q * q + 2) / 6,
        (q * q + 2) / (2 * q * q),
        (q * q - 2 * q + 2) / (q * q * (q - 1)),
        1 / (q * q),
    ]


def double_endpoint_coefficient(split: tuple[int, ...]) -> float:
    energies = endpoint_energies()
    left = split[0]
    right = split[3]
    return min(
        sqrt(energies[left] * energies[right]),
        sqrt(energies[3 - left] * energies[3 - right]),
    )


def adjacent_degrees() -> tuple[list[int], list[int], float]:
    q = Q
    endpoint = [
        q * comb(q, 3) + q * q * (q - 1) * comb(q, 2),
        comb(q - 1, 2) + (q - 1) * comb(q, 2) + q * (q - 1) ** 2,
        q * q - 2,
        1,
    ]
    middle = [q * q * (q - 1) ** 2, 3 * (q - 1) ** 2, 2 * (q - 1), 1]
    maximum_squared = (q + 2) ** 2 / (
        q * q * (q - 1) ** 2 * (q - 2) ** 2
    )
    return endpoint, middle, maximum_squared


def adjacent_coefficient(split: tuple[int, ...]) -> float:
    endpoint, middle, maximum_squared = adjacent_degrees()
    left = split[0]
    right = split[1]
    return min(
        sqrt(maximum_squared * endpoint[left] * middle[right]) / (Q - 1),
        sqrt(
            maximum_squared * endpoint[3 - left] * middle[3 - right]
        )
        / (Q * (Q - 1)),
    )


def deterministic_cut_sum(
    profile: tuple[int, ...],
    occupation: tuple[int, ...],
    coefficient,
) -> float:
    result = 0.0
    for split in product(*(range(degree + 1) for degree in profile)):
        ket = 1
        bra = 1
        for block in range(4):
            ket *= comb(occupation[block], split[block])
            bra *= comb(
                occupation[block], profile[block] - split[block]
            )
        result += coefficient(split) * sqrt(ket * bra)
    return result


def main() -> None:
    double_endpoint = deterministic_cut_sum(
        (3, 1, 1, 3),
        (2, 1, 1, 2),
        double_endpoint_coefficient,
    )
    adjacent = deterministic_cut_sum(
        (3, 3, 1, 1),
        (2, 2, 1, 1),
        adjacent_coefficient,
    )
    attenuated_endpoint = BETA**8 * double_endpoint
    attenuated_adjacent = BETA**8 * adjacent

    minimal_coefficient = 2337 / 256 + 3 * sqrt(2) / 8
    promise_loss = 2 * 1288991 / 94064415
    committed = BETA**4 * minimal_coefficient / Q + promise_loss
    remaining_margin = 1 / 3 - committed
    required_endpoint_factor = remaining_margin / attenuated_endpoint

    if not abs(double_endpoint - 2.456327355938791) < 2e-13:
        raise AssertionError(("double-endpoint cut sum", double_endpoint))
    if not abs(adjacent - 0.10800960164732769) < 2e-13:
        raise AssertionError(("adjacent cut sum", adjacent))
    if required_endpoint_factor >= 0.281:
        raise AssertionError(("required improvement factor", required_endpoint_factor))

    print(
        "fixed-split occupation barrier confirmed: "
        f"double_endpoint_raw={double_endpoint:.12g},"
        f"double_endpoint_beta8={attenuated_endpoint:.12g},"
        f"adjacent_raw={adjacent:.12g},"
        f"adjacent_beta8={attenuated_adjacent:.12g},"
        f"minimal_plus_promise={committed:.12g},"
        f"remaining_TV_margin={remaining_margin:.12g},"
        f"required_endpoint_factor_below={required_endpoint_factor:.12g}"
    )


if __name__ == "__main__":
    main()
