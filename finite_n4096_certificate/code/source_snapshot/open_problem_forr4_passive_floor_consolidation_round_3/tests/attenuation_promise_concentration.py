#!/usr/bin/env python3
"""Regression for the Kearns--Saul attenuation promise bound."""

from __future__ import annotations

from math import exp, log, sqrt
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from attenuation_promise_concentration import (
    biased_sign_lower_log_mgf,
    biased_rademacher_proxy,
    euclidean_promise_concentration,
    extended_euclidean_promise_concentration,
    hybrid_euclidean_promise_concentration,
    promise_concentration,
    two_split_euclidean_gate,
)
from promise_tail_monte_carlo import exact_plant, sylvester_sign


def scalar_mgf_check(beta: float) -> None:
    proxy = biased_rademacher_proxy(beta)
    positive = (1 + beta) / 2
    for value in np.linspace(-18, 18, 1441):
        mgf = positive * exp(value * (1 - beta)) + (
            1 - positive
        ) * exp(value * (-1 - beta))
        upper = exp(proxy * value * value / 2)
        if mgf > upper * (1 + 2e-13):
            raise AssertionError(("Kearns--Saul mgf", beta, value, mgf, upper))


def monotone_branch_check(beta: float) -> None:
    """Regression for the exact branch used by Euclidean packing."""

    kearns_saul_tilt = log((1 + beta) / (1 - beta))
    values = np.linspace(kearns_saul_tilt / 2000, kearns_saul_tilt, 2000)
    ratios = np.asarray(
        [biased_sign_lower_log_mgf(beta, value) / value**2 for value in values]
    )
    if np.min(np.diff(ratios)) < -3e-11:
        raise AssertionError(("nonmonotone exact MGF branch", beta))
    endpoint = biased_sign_lower_log_mgf(beta, kearns_saul_tilt)
    if not np.isclose(endpoint, beta * kearns_saul_tilt, atol=2e-14):
        raise AssertionError(("Kearns--Saul equality tilt", beta, endpoint))


def euclidean_packing_check(beta: float) -> None:
    """Compare the packed bound with exact products for signed weights."""

    rng = np.random.default_rng(19073)
    positive = (1 + beta) / 2
    negative = (1 - beta) / 2
    kearns_saul_tilt = log((1 + beta) / (1 - beta))
    for width in (1, 2, 5, 17):
        for _ in range(30):
            coefficients = rng.normal(size=width)
            norm = float(np.linalg.norm(coefficients))
            tilt = rng.uniform(0.02, 1) * kearns_saul_tilt / norm
            exact = 0.0
            for coefficient in coefficients:
                terms = (
                    log(positive) - tilt * coefficient * (1 - beta),
                    log(negative) + tilt * coefficient * (1 + beta),
                )
                exact += float(np.logaddexp(*terms))
            packed = biased_sign_lower_log_mgf(beta, tilt * norm)
            if exact > packed + 2e-13:
                raise AssertionError(
                    ("Euclidean MGF packing", beta, width, exact, packed)
                )


def extended_pair_packing_check(beta: float) -> None:
    """Protect the scalar inequality behind the two-split extension."""

    branch = log((1 + beta) / (1 - beta))
    gate = two_split_euclidean_gate(beta)
    if not branch < gate < sqrt(2) * branch:
        raise AssertionError(("two-split gate", beta, branch, gate))
    if not np.isclose(
        biased_sign_lower_log_mgf(beta, gate),
        2 * biased_sign_lower_log_mgf(beta, gate / sqrt(2)),
        atol=2e-14,
    ):
        raise AssertionError(("two-split equality", beta, gate))
    for total in np.linspace(branch, gate, 31):
        packed = biased_sign_lower_log_mgf(beta, total)
        for first in np.linspace(total / sqrt(2), total, 121):
            second = sqrt(max(0.0, total**2 - first**2))
            pair = (
                biased_sign_lower_log_mgf(beta, first)
                + biased_sign_lower_log_mgf(beta, second)
            )
            if pair > packed + 3e-13:
                raise AssertionError(
                    (
                        "extended pair packing",
                        beta,
                        total,
                        first,
                        pair,
                        packed,
                    )
                )
        for count in range(2, 33):
            equal = count * biased_sign_lower_log_mgf(
                beta,
                total / sqrt(count),
            )
            if equal > packed + 3e-13:
                raise AssertionError(
                    ("extended equal packing", beta, total, count, equal, packed)
                )


def reverse_martingale_norm_check(beta: float) -> None:
    """Check the four exact conditional coefficient norms at q=4."""

    rng = np.random.default_rng(28121)
    order = 4
    dimension = order**2
    plant = exact_plant(sylvester_sign(order), rng)
    hadamard = sylvester_sign(dimension) / sqrt(dimension)
    chain = np.einsum(
        "ij,jk,kl->ijkl", hadamard, hadamard, hadamard
    ) / dimension
    for _ in range(8):
        noise = tuple(rng.choice((-1, 1), dimension) for _ in range(4))
        for level in range(4):
            weights = tuple(
                plant[block]
                * (beta if block < level else noise[block])
                for block in range(4)
            )
            if level == 0:
                coefficient = plant[0] * np.einsum(
                    "ijkl,j,k,l->i", chain, weights[1], weights[2], weights[3]
                )
            elif level == 1:
                coefficient = plant[1] * np.einsum(
                    "ijkl,i,k,l->j", chain, weights[0], weights[2], weights[3]
                )
            elif level == 2:
                coefficient = plant[2] * np.einsum(
                    "ijkl,i,j,l->k", chain, weights[0], weights[1], weights[3]
                )
            else:
                coefficient = plant[3] * np.einsum(
                    "ijkl,i,j,k->l", chain, weights[0], weights[1], weights[2]
                )
            expected = beta**level / sqrt(dimension)
            if not np.isclose(np.linalg.norm(coefficient), expected, atol=2e-14):
                raise AssertionError(
                    (
                        "reverse martingale coefficient norm",
                        level,
                        np.linalg.norm(coefficient),
                        expected,
                    )
                )


def main() -> None:
    for beta in (0.1, 0.5, 3 / 4, 313 / 400, 5 / 6, 0.95):
        scalar_mgf_check(beta)
        monotone_branch_check(beta)
        euclidean_packing_check(beta)
    for beta in (3 / 4, 31 / 40, 0.77985, 79 / 100):
        extended_pair_packing_check(beta)
    reverse_martingale_norm_check(0.780899845855353)

    result = promise_concentration(1024, 313 / 400)
    if not 0.0157 < result.two_hypothesis_loss < 0.0158:
        raise AssertionError(result)
    if not result.mean > 1 / 4:
        raise AssertionError(result)
    repaired = euclidean_promise_concentration(1024, 0.780899845855353)
    if not 0.01818 < repaired.two_hypothesis_loss < 0.01820:
        raise AssertionError(repaired)
    if not np.isclose(
        repaired.chernoff_tilt,
        repaired.kearns_saul_tilt * sqrt(repaired.dimension),
        atol=1e-14,
    ):
        raise AssertionError(repaired)
    hybrid = hybrid_euclidean_promise_concentration(
        1024,
        0.780899845855353,
    )
    if hybrid.globalized_blocks != 1:
        raise AssertionError(("hybrid globalized blocks", hybrid))
    if not (
        hybrid.kearns_saul_tilt
        < hybrid.scalar_chernoff_tilt
        < hybrid.kearns_saul_tilt / hybrid.beta
    ):
        raise AssertionError(("hybrid tilt branch", hybrid))
    if not np.isclose(
        hybrid.two_hypothesis_loss,
        0.0175468358426,
        atol=3e-12,
    ):
        raise AssertionError(("hybrid promise loss", hybrid))
    if hybrid.two_hypothesis_loss >= repaired.two_hypothesis_loss:
        raise AssertionError(("hybrid promise improvement", hybrid, repaired))
    extended = extended_euclidean_promise_concentration(
        1024,
        0.779844397452418,
    )
    same_beta_hybrid = hybrid_euclidean_promise_concentration(
        1024,
        extended.beta,
    )
    if not np.isclose(
        extended.two_hypothesis_loss,
        0.0200666604394,
        atol=3e-12,
    ):
        raise AssertionError(("extended promise loss", extended))
    if not (
        extended.kearns_saul_tilt
        < extended.scalar_chernoff_tilt
        < extended.two_split_gate
    ):
        raise AssertionError(("extended promise branch", extended))
    if extended.two_hypothesis_loss >= same_beta_hybrid.two_hypothesis_loss:
        raise AssertionError(
            ("extended promise improvement", extended, same_beta_hybrid)
        )
    print(
        "attenuation promise concentration passed: "
        f"beta={result.beta:.12g},mean={result.mean:.12g},"
        f"sign_proxy={result.rademacher_proxy:.12g},"
        f"chain_proxy={result.four_chain_proxy:.12g},"
        f"conditioning_loss={result.two_hypothesis_loss:.12g},"
        f"repaired_loss={repaired.two_hypothesis_loss:.12g},"
        f"hybrid_loss={hybrid.two_hypothesis_loss:.12g},"
        f"extended_loss={extended.two_hypothesis_loss:.12g},"
        f"same_beta_hybrid_loss={same_beta_hybrid.two_hypothesis_loss:.12g}"
    )


if __name__ == "__main__":
    main()
