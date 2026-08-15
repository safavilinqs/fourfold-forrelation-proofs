#!/usr/bin/env python3
"""Promise bounds for the independently attenuated exact plant."""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, log, sqrt

from scipy.optimize import brentq, minimize_scalar


@dataclass(frozen=True)
class PromiseConcentration:
    dimension: int
    beta: float
    mean: float
    promise_gap: float
    rademacher_proxy: float
    four_chain_proxy: float
    one_sided_failure: float
    two_hypothesis_loss: float


@dataclass(frozen=True)
class EuclideanPromiseConcentration:
    """Finite-tilt four-block concentration certificate.

    The certificate uses the monotonic branch of the exact biased-sign log
    MGF rather than replacing every martingale difference by its global
    Kearns--Saul parabola.
    """

    dimension: int
    beta: float
    mean: float
    promise_gap: float
    kearns_saul_tilt: float
    chernoff_tilt: float
    log_one_sided_failure: float
    one_sided_failure: float
    two_hypothesis_loss: float


@dataclass(frozen=True)
class HybridEuclideanPromiseConcentration:
    """Optimized exact/global packing across the four martingale blocks."""

    dimension: int
    beta: float
    mean: float
    promise_gap: float
    kearns_saul_tilt: float
    scalar_chernoff_tilt: float
    chernoff_tilt: float
    globalized_blocks: int
    log_one_sided_failure: float
    one_sided_failure: float
    two_hypothesis_loss: float


@dataclass(frozen=True)
class ExtendedEuclideanPromiseConcentration:
    dimension: int
    beta: float
    mean: float
    promise_gap: float
    kearns_saul_tilt: float
    two_split_gate: float
    scalar_chernoff_tilt: float
    chernoff_tilt: float
    log_one_sided_failure: float
    one_sided_failure: float
    two_hypothesis_loss: float


def biased_rademacher_proxy(beta: float) -> float:
    """Optimal Kearns--Saul subgaussian proxy for a biased sign."""

    if not 0 < beta < 1:
        raise ValueError(beta)
    return 2 * beta / log((1 + beta) / (1 - beta))


def biased_sign_lower_log_mgf(beta: float, value: float) -> float:
    """Return log E exp[-value (eta-beta)] for E eta=beta.

    ``value`` is restricted to the lower-tail direction.  The stable
    two-term log-sum-exp form avoids overflow at large tilts.
    """

    if not 0 < beta < 1:
        raise ValueError(beta)
    if value < 0:
        raise ValueError(value)
    positive = (1 + beta) / 2
    negative = (1 - beta) / 2
    terms = (
        log(positive) - value * (1 - beta),
        log(negative) + value * (1 + beta),
    )
    pivot = max(terms)
    return pivot + log(sum(exp(term - pivot) for term in terms))


def promise_concentration(
    dimension: int, beta: float, promise_boundary: float = 1 / 4
) -> PromiseConcentration:
    if dimension <= 0:
        raise ValueError(dimension)
    mean = beta**4
    gap = mean - promise_boundary
    if gap <= 0:
        raise ValueError(("mean does not clear the promise", mean))
    sign_proxy = biased_rademacher_proxy(beta)
    # Two sequential bilinear Hadamard contractions.  Conditional endpoint
    # noise has proxy kappa(1+beta^2)/N; the middle two-block form has the
    # same proxy and is multiplied by beta^2.
    chain_proxy = (
        sign_proxy * (1 + beta**2) * (1 + beta**4) / dimension
    )
    failure = exp(-(gap**2) / (2 * chain_proxy))
    return PromiseConcentration(
        dimension=dimension,
        beta=beta,
        mean=mean,
        promise_gap=gap,
        rademacher_proxy=sign_proxy,
        four_chain_proxy=chain_proxy,
        one_sided_failure=failure,
        two_hypothesis_loss=2 * failure,
    )


def euclidean_promise_concentration(
    dimension: int, beta: float, promise_boundary: float = 1 / 4
) -> EuclideanPromiseConcentration:
    """Sharpen the promise tail using exact finite-tilt Euclidean packing.

    For a centered sign of mean ``beta``, Schlemm's strict-unimodality
    theorem for the Bernoulli log-MGF implies that

        psi(s) / s**2 is increasing on 0 <= s <= L,
        L = log((1+beta)/(1-beta)).

    Consequently, a weighted centered-sign sum with coefficient l2 norm A
    has log-MGF at most psi(t A) whenever t A <= L.  The four reverse Doob
    differences of the planted chain have conditional coefficient norms
    beta**j / sqrt(dimension), j=0,1,2,3.  Iterating the exact bound and
    choosing the explicit valid tilt t=L sqrt(dimension) gives the returned
    one-sided certificate.
    """

    if dimension <= 0:
        raise ValueError(dimension)
    if not 0 < beta < 1:
        raise ValueError(beta)
    mean = beta**4
    gap = mean - promise_boundary
    if gap <= 0:
        raise ValueError(("mean does not clear the promise", mean))
    scalar_tilt = log((1 + beta) / (1 - beta))
    chernoff_tilt = scalar_tilt * sqrt(dimension)
    log_failure = -chernoff_tilt * gap + sum(
        biased_sign_lower_log_mgf(beta, scalar_tilt * beta**level)
        for level in range(4)
    )
    failure = min(1.0, exp(log_failure))
    return EuclideanPromiseConcentration(
        dimension=dimension,
        beta=beta,
        mean=mean,
        promise_gap=gap,
        kearns_saul_tilt=scalar_tilt,
        chernoff_tilt=chernoff_tilt,
        log_one_sided_failure=log_failure,
        one_sided_failure=failure,
        two_hypothesis_loss=2 * failure,
    )


def hybrid_euclidean_promise_concentration(
    dimension: int,
    beta: float,
    promise_boundary: float = 1 / 4,
) -> HybridEuclideanPromiseConcentration:
    """Optimize finite-branch and global packing block by block.

    Write s=t/sqrt(dimension).  Martingale block j has Euclidean
    coefficient norm beta**j/sqrt(dimension).  If s*beta**j <= L the
    finite-tilt theorem gives the exact packed term psi(s*beta**j).
    Beyond that branch, the global Kearns--Saul parabola gives
    kappa*(s*beta**j)**2/2.  The two formulas agree with matching first
    derivative at L.  Optimizing over every branch interval is therefore
    a valid Chernoff certificate, not an extrapolation of the finite-tilt
    inequality.
    """

    if dimension <= 0:
        raise ValueError(dimension)
    if not 0 < beta < 1:
        raise ValueError(beta)
    mean = beta**4
    gap = mean - promise_boundary
    if gap <= 0:
        raise ValueError(("mean does not clear the promise", mean))
    branch_tilt = log((1 + beta) / (1 - beta))
    proxy = biased_rademacher_proxy(beta)

    def log_failure(scalar_tilt: float) -> float:
        result = -scalar_tilt * sqrt(dimension) * gap
        for level in range(4):
            value = scalar_tilt * beta**level
            if value <= branch_tilt:
                result += biased_sign_lower_log_mgf(beta, value)
            else:
                result += proxy * value**2 / 2
        return result

    boundaries = [0.0] + [
        branch_tilt / beta**level
        for level in range(4)
    ]
    # The all-global objective is quadratic and has unconstrained minimizer
    # shown below.  Include it as the finite right endpoint of the search.
    global_variance = proxy * sum(beta ** (2 * level) for level in range(4))
    global_minimizer = sqrt(dimension) * gap / global_variance
    boundaries.append(max(boundaries[-1], global_minimizer))
    candidates = set(boundaries)
    for left, right in zip(boundaries[:-1], boundaries[1:], strict=True):
        if right <= left:
            continue
        optimum = minimize_scalar(
            log_failure,
            bounds=(left, right),
            method="bounded",
            options={"xatol": 1e-13},
        )
        candidates.add(float(optimum.x))
    scalar_tilt = min(candidates, key=log_failure)
    log_bound = log_failure(scalar_tilt)
    failure = min(1.0, exp(log_bound))
    globalized = sum(
        scalar_tilt * beta**level > branch_tilt
        for level in range(4)
    )
    return HybridEuclideanPromiseConcentration(
        dimension=dimension,
        beta=beta,
        mean=mean,
        promise_gap=gap,
        kearns_saul_tilt=branch_tilt,
        scalar_chernoff_tilt=scalar_tilt,
        chernoff_tilt=scalar_tilt * sqrt(dimension),
        globalized_blocks=globalized,
        log_one_sided_failure=log_bound,
        one_sided_failure=failure,
        two_hypothesis_loss=2 * failure,
    )


def two_split_euclidean_gate(beta: float) -> float:
    """Return the first tilt where two equal coefficients beat one.

    For the lower-tail biased-sign log-MGF `psi`, fixed-Euclidean-norm
    extrema have equal nonzero coefficient magnitudes.  Up to this gate,
    one coefficient is the worst case and the packed block costs
    `psi(s)` even though `s` is slightly beyond the Kearns--Saul
    equality tilt.
    """

    if not 0 < beta < 1:
        raise ValueError(beta)
    branch_tilt = log((1 + beta) / (1 - beta))

    def difference(value: float) -> float:
        return (
            2
            * biased_sign_lower_log_mgf(
                beta,
                value / sqrt(2),
            )
            - biased_sign_lower_log_mgf(beta, value)
        )

    upper = sqrt(2) * branch_tilt
    if difference(branch_tilt) >= 0 or difference(upper) <= 0:
        raise AssertionError(
            (
                "two-split gate bracket",
                beta,
                difference(branch_tilt),
                difference(upper),
            )
        )
    return float(brentq(difference, branch_tilt, upper, xtol=1e-14))


def extended_euclidean_promise_concentration(
    dimension: int,
    beta: float,
    promise_boundary: float = 1 / 4,
) -> ExtendedEuclideanPromiseConcentration:
    """Optimize exact Euclidean block packing through the two-split gate."""

    if dimension <= 0:
        raise ValueError(dimension)
    if not 0 < beta < 1:
        raise ValueError(beta)
    mean = beta**4
    gap = mean - promise_boundary
    if gap <= 0:
        raise ValueError(("mean does not clear the promise", mean))
    branch_tilt = log((1 + beta) / (1 - beta))
    gate = two_split_euclidean_gate(beta)

    def log_failure(scalar_tilt: float) -> float:
        return (
            -scalar_tilt * sqrt(dimension) * gap
            + sum(
                biased_sign_lower_log_mgf(
                    beta,
                    scalar_tilt * beta**level,
                )
                for level in range(4)
            )
        )

    optimum = minimize_scalar(
        log_failure,
        bounds=(0, gate),
        method="bounded",
        options={"xatol": 1e-13},
    )
    candidates = (0.0, float(optimum.x), gate)
    scalar_tilt = min(candidates, key=log_failure)
    log_bound = log_failure(scalar_tilt)
    failure = min(1.0, exp(log_bound))
    return ExtendedEuclideanPromiseConcentration(
        dimension=dimension,
        beta=beta,
        mean=mean,
        promise_gap=gap,
        kearns_saul_tilt=branch_tilt,
        two_split_gate=gate,
        scalar_chernoff_tilt=scalar_tilt,
        chernoff_tilt=scalar_tilt * sqrt(dimension),
        log_one_sided_failure=log_bound,
        one_sided_failure=failure,
        two_hypothesis_loss=2 * failure,
    )


def main() -> None:
    for beta in (3 / 4, 31 / 40, 313 / 400, 79 / 100, 5 / 6):
        result = promise_concentration(1024, beta)
        print(
            f"beta={beta:.12g},mean={result.mean:.12g},"
            f"proxy={result.four_chain_proxy:.12g},"
            f"one_sided_failure={result.one_sided_failure:.12g},"
            f"conditioning_loss={result.two_hypothesis_loss:.12g}"
        )
        sharpened = euclidean_promise_concentration(1024, beta)
        print(
            f"beta={beta:.12g},euclidean_tilt={sharpened.chernoff_tilt:.12g},"
            f"euclidean_one_sided_failure="
            f"{sharpened.one_sided_failure:.12g},"
            f"euclidean_conditioning_loss="
            f"{sharpened.two_hypothesis_loss:.12g}"
        )
        hybrid = hybrid_euclidean_promise_concentration(1024, beta)
        print(
            f"beta={beta:.12g},hybrid_scalar_tilt="
            f"{hybrid.scalar_chernoff_tilt:.12g},"
            f"globalized_blocks={hybrid.globalized_blocks},"
            f"hybrid_conditioning_loss={hybrid.two_hypothesis_loss:.12g}"
        )
        extended = extended_euclidean_promise_concentration(1024, beta)
        print(
            f"beta={beta:.12g},extended_gate="
            f"{extended.two_split_gate:.12g},"
            f"extended_scalar_tilt={extended.scalar_chernoff_tilt:.12g},"
            f"extended_conditioning_loss="
            f"{extended.two_hypothesis_loss:.12g}"
        )


if __name__ == "__main__":
    main()
