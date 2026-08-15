#!/usr/bin/env python3
"""Translation-twirled arbitrary-diagonal bounds for one endpoint cubic.

The profile is (3,1,1,1), or its reversal.  For an internal 1|2 split of
the cubic support, translation twirling leaves only the endpoint-pair XOR
law.  Exact endpoint weight-matrix spectra then give closed coefficients for
balanced and extreme placements of the three singleton chain coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext


@dataclass(frozen=True)
class EndpointCubicBound:
    order: int
    extreme_singletons: Decimal
    balanced_singletons: Decimal


def bound(order: int) -> EndpointCubicBound:
    if order < 2 or order & (order - 1):
        raise ValueError(("order must be a power of two", order))
    with localcontext() as context:
        context.prec = 70
        q = Decimal(order)
        n = q - 1
        large = q * q - 2 * q + 2
        common_root = (n * (large * large + q**3)).sqrt()

        # With all three singleton coordinates on one side, the endpoint
        # pair is the only variable on the other side.  The two endpoint
        # difference classes give this exact nuclear sum.  The reverse
        # extreme orientation has a smaller rank--Frobenius bound.
        extreme = 2 * common_root / q**5

        # With one or two singleton coordinates beside the endpoint
        # singleton, the pair-difference column Grams are orthogonal.  Their
        # per-type nuclear sums are q*sqrt(2)*A and q^2*sqrt(2), producing an
        # extra factor q relative to the extreme placement.
        balanced = 2 * common_root / q**4
    return EndpointCubicBound(
        order=order,
        extreme_singletons=extreme,
        balanced_singletons=balanced,
    )


def main() -> None:
    for order in (2, 4, 8, 16, 32):
        result = bound(order)
        print(
            f"q={order},N={order * order},"
            f"extreme={result.extreme_singletons},"
            f"balanced={result.balanced_singletons}"
        )


if __name__ == "__main__":
    main()
