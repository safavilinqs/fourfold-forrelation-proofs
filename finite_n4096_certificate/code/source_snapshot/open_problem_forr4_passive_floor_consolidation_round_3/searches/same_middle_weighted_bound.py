#!/usr/bin/env python3
"""Exact arbitrary-diagonal bounds for same-middle endpoint placements.

For the profile (3,1,1,3), put both middle singleton occurrences on the
same side of the ket/bra cut.  Translation twirling leaves only the XOR
differences of the two endpoint pairs.  The resulting weighted matrices
have an explicit Gram reduction, so the two endpoint orientations can be
optimized exactly for every power-of-two order q.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext


@dataclass(frozen=True)
class SameMiddleBound:
    order: int
    endpoint_orbit_large: int
    equal_endpoint_orientation: Decimal
    mixed_endpoint_orientation: Decimal
    hybrid_pair_column: Decimal
    hybrid_whole_column: Decimal
    hybrid_upper: Decimal


def bound(order: int) -> SameMiddleBound:
    """Return the exact coefficients for equal and mixed orientations."""

    if order < 2 or order & (order - 1):
        raise ValueError(("order must be a power of two", order))
    with localcontext() as context:
        context.prec = 70
        q = Decimal(order)
        n = q - 1
        orbit_large = order * order - 2 * order + 2

        # S_V=A/(2(q-1)) and S_perp=1/(q-1) are the squared
        # endpoint-pair column energies at fixed XOR difference.
        equal = Decimal(orbit_large) / (n * q**3)

        # The endpoint weight matrix W_x has nuclear norm sqrt(2) A for a
        # vertical difference and sqrt(2) q otherwise.  Optimizing the two
        # remaining difference laws gives the closed expression below.
        mixed = (
            2
            * (
                Decimal(orbit_large)
                * (Decimal(orbit_large) ** 2 + q**3)
            ).sqrt()
            / q**5
        )
        # Exactly one endpoint is split 1|2 and the other remains a whole
        # cubic block.  If the pair is the only column variable, the endpoint
        # nuclear spectrum gives the first value.  In the transposed
        # placement, rank--Frobenius on the whole cubic Bessel frame gives
        # the second.  The latter is the safe orientation-independent bound.
        hybrid_pair_column = (
            2
            * (
                n
                * (Decimal(orbit_large) ** 2 + q**3)
            ).sqrt()
            / q**5
        )
        hybrid_whole_column = (
            (Decimal(orbit_large) / n).sqrt() / q**2
        )
        hybrid_upper = max(hybrid_pair_column, hybrid_whole_column)
    return SameMiddleBound(
        order=order,
        endpoint_orbit_large=orbit_large,
        equal_endpoint_orientation=equal,
        mixed_endpoint_orientation=mixed,
        hybrid_pair_column=hybrid_pair_column,
        hybrid_whole_column=hybrid_whole_column,
        hybrid_upper=hybrid_upper,
    )


def deterministic_ledger(
    alternating_equal: Decimal = Decimal("0.010905"),
    alternating_mixed: Decimal = Decimal("0.0306880103312"),
) -> tuple[Decimal, Decimal, Decimal]:
    """Safe q=32 ledger for the deterministic occupation (2,1,1,2)."""

    result = bound(32)
    with localcontext() as context:
        context.prec = 70
        # There are two equal endpoint orientations and two mixed endpoint
        # orientations.  Their combined occurrence mass is four in either
        # class, and each class has two same-middle and two alternating
        # middle placements.
        raw = 4 * (
            2 * result.equal_endpoint_orientation
            + 2 * alternating_equal
        ) + 4 * (
            2 * result.mixed_endpoint_orientation
            + 2 * alternating_mixed
        )
        attenuated = (Decimal(5) / 6) ** 8 * raw
        remaining_margin = Decimal("0.160358131958")
    return raw, attenuated, remaining_margin


def main() -> None:
    for order in (2, 4, 8, 16, 32):
        result = bound(order)
        print(
            f"q={order},N={order * order},"
            f"equal_orientation={result.equal_endpoint_orientation},"
            f"mixed_orientation={result.mixed_endpoint_orientation},"
            f"hybrid_pair_column={result.hybrid_pair_column},"
            f"hybrid_whole_column={result.hybrid_whole_column}"
        )
    raw, attenuated, margin = deterministic_ledger()
    print(
        "same-middle-refined double-endpoint ledger: "
        f"raw={raw},attenuated={attenuated},"
        f"available_margin={margin},slack={margin-attenuated}"
    )


if __name__ == "__main__":
    main()
