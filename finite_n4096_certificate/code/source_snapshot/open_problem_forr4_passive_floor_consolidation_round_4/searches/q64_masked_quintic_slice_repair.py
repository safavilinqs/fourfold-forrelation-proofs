#!/usr/bin/env python3
"""Repair 54 masked all-record-one quintic entries by physical factorizations.

The physical disjointness mask is retained.  For a fixed occurrence row,
every surviving column completes its partial supports.  Exact degree-five
singleton slice energies bound the squared norm of that complete physical
row.  Factoring the kernel through its row vectors then gives an
arbitrary-diagonal trace-norm coefficient equal to the square root of the
better row/complement energy.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from json import dumps
from math import comb, inf, nextafter, sqrt
from pathlib import Path

from q64_noncubic_recovered_universal_insertion import (
    recovered_universal_entries,
    universal_noncubic_entries,
)
from q64_universal_double_cubic_insertion import double_cubic_entries
from q64_universal_multicubic_insertion import multicubic_entries
from q64_universal_septimic_insertion import septimic_entries


ROOT = Path(__file__).resolve().parents[1]
ORDER = 64
PROFILE_FAMILIES = (
    (1, 5, 1, 5),
    (5, 1, 1, 5),
    (5, 1, 5, 1),
)

Profile = tuple[int, ...]
Split = tuple[int, ...]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class MaskedQuinticSliceRepair:
    order: int
    dimension: int
    coefficient_one_dependent_entries: int
    candidate_entries: int
    repaired_entries: int
    remaining_quarantined_entries: int
    repaired_profiles: int
    bi_record_quintic_tail_bound: float
    distinct_squared_coefficients: tuple[float, ...]
    minimum_coefficient: float
    maximum_coefficient: float
    maximum_squared_coefficient: float


def coefficient_one_dependent_entries() -> frozenset[ProfileSplit]:
    classes = (
        septimic_entries(),
        multicubic_entries(),
        double_cubic_entries(),
        universal_noncubic_entries(),
        recovered_universal_entries(),
    )
    result = frozenset(entry for entries in classes for entry in entries)
    if len(result) != sum(map(len, classes)):
        raise AssertionError("coefficient-one inventories overlap")
    return result


def exact_quintic_singleton_slice_energies(
    order: int,
) -> tuple[Fraction, ...]:
    """Exact rational form of the accepted degree-five slice identities."""

    q = order
    if q < 4 or q & (q - 1):
        raise ValueError(("quintic slices require power-of-two order", q))
    w0 = Fraction(1, q**2)
    w1 = Fraction(1, q**2 * (q - 1) ** 2)
    w2 = Fraction(4, q**2 * (q - 1) ** 2 * (q - 2) ** 2)

    zero_xor_four_sets = q * (q - 1) * (q - 2) // 24
    count_5 = q * comb(q, 5)
    count_41 = q**2 * (q - 1) * comb(q, 4)
    high_41 = q**2 * (q - 1) * zero_xor_four_sets
    count_32 = q * (q - 1) * comb(q, 3) * comb(q, 2)
    count_221 = q**2 * comb(q - 1, 2) * comb(q, 2) ** 2
    high_221 = q**2 * comb(q - 1, 2) * (q - 1) * q**2 // 4
    slice_0 = (
        (count_5 + high_41) * w0
        + (count_41 - high_41 + count_32 + high_221) * w1
        + (count_221 - high_221) * w2
    )
    slice_1 = Fraction(5, q**2) * slice_0

    half_nonzero_xors = q // 2 - 1
    slice_2 = (
        comb(q - 2, 3) * w0
        + q
        * (q - 1)
        * (
            half_nonzero_xors * w0
            + (comb(q - 2, 2) - half_nonzero_xors) * w1
        )
        + (
            (q - 2) * (q - 1) * comb(q, 2)
            + (q - 1) * comb(q, 3)
        )
        * w1
        + (q - 1)
        * (q - 2)
        * q
        * (Fraction(q, 2) * w1 + (comb(q, 2) - q // 2) * w2)
    )
    slice_3 = (
        comb(q - 3, 2) * w0
        + q * (q - 1) * (w0 + (q - 4) * w1)
        + (q - 1) * comb(q, 2) * w1
    )
    return (
        slice_0,
        slice_1,
        slice_2,
        slice_3,
        Fraction(q**2 - 4, q**2),
        Fraction(1, q**2),
    )


def complement(profile: Profile, split: Split) -> Split:
    return tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )


def row_energy_bound(profile: Profile, split: Split) -> Fraction:
    """Bound one full physical row's squared energy, masks included."""

    q = ORDER
    energy = exact_quintic_singleton_slice_energies(q)
    if profile == (5, 1, 1, 5):
        # The central singleton--singleton link has squared magnitude q^-2.
        # A variable singleton endpoint contributes q^2 copies of its exact
        # quintic slice; a fixed singleton contributes one copy.
        exponent = 2 * (1 - split[1] - split[2])
        scale = (
            Fraction(q**exponent, 1)
            if exponent >= 0
            else Fraction(1, q ** (-exponent))
        )
        return scale * energy[split[0]] * energy[split[3]]
    if profile == (5, 1, 5, 1):
        return row_energy_bound((1, 5, 1, 5), tuple(reversed(split)))
    if profile != (1, 5, 1, 5):
        raise ValueError(("not a masked quintic slice profile", profile))

    fixed_singletons = split[0] + split[2]
    if fixed_singletons == 0:
        # Every internal quintic has record one on both adjacent axes.
        # This excludes the degree-five exceptional supports with singleton
        # tail energy one, leaving tail energy at most (q-1)^-2.
        middle_scale = Fraction(q**2, (q - 1) ** 2)
    elif fixed_singletons == 1:
        middle_scale = Fraction(1, (q - 1) ** 2)
    else:
        middle_scale = Fraction(1, q**2)
    return energy[split[1]] * energy[split[3]] * middle_scale


def squared_coefficient(entry: ProfileSplit) -> Fraction:
    profile, split = entry
    bounds = [
        row_energy_bound(profile, split),
        row_energy_bound(profile, complement(profile, split)),
    ]
    if profile == (5, 1, 1, 5):
        bounds.append(separated_tensor_squared_bound(split))
        if split[1] != split[2]:
            bounds.append(quintic_hadamard_chain_squared_bound(split))
    return min(bounds)


def endpoint_link_squared_factor(selected: int, singleton_selected: int) -> Fraction:
    """Best row/column feature-norm product for one masked endpoint link."""

    q = ORDER
    energy = exact_quintic_singleton_slice_energies(q)
    row_energy = q ** (2 * (1 - singleton_selected)) * energy[selected]
    column_energy = q ** (2 * singleton_selected) * energy[5 - selected]
    return min(row_energy, column_energy)


def separated_tensor_squared_bound(split: Split) -> Fraction:
    """Tensor-factor bound for the profile (5,1,1,5)."""

    q = ORDER
    left = endpoint_link_squared_factor(split[0], split[1])
    right = endpoint_link_squared_factor(split[3], split[2])
    central = Fraction(1, q**2) if split[1] == split[2] else Fraction(1, 1)
    return left * right * central


def cross_distinctness_factor_upper(selected: int) -> int:
    """Integer upper bound for the selected/complement quintic mask."""

    depth = min(selected, 5 - selected)
    if depth == 0:
        return 1
    if depth == 1:
        return 3
    if depth == 2:
        return 6
    raise ValueError(("not a quintic split", selected))


def quintic_hadamard_chain_squared_bound(split: Split) -> Fraction:
    """Masked completion bound retaining the central Walsh collapse."""

    if split[1] == split[2]:
        raise ValueError(("central singletons must be opposite", split))
    factor = (
        cross_distinctness_factor_upper(split[0])
        * cross_distinctness_factor_upper(split[3])
    )
    return Fraction(factor, ORDER) ** 2


def candidate_entries() -> tuple[ProfileSplit, ...]:
    affected = coefficient_one_dependent_entries()
    return tuple(
        sorted(entry for entry in affected if entry[0] in PROFILE_FAMILIES)
    )


def repaired_entries() -> tuple[ProfileSplit, ...]:
    return tuple(
        entry
        for entry in candidate_entries()
        if squared_coefficient(entry) <= 1
    )


def outward_coefficient(entry: ProfileSplit) -> float:
    exact = squared_coefficient(entry)
    result = sqrt(float(exact))
    while Fraction.from_float(result) ** 2 < exact:
        result = nextafter(result, inf)
    return result


def coefficient_map() -> dict[ProfileSplit, float]:
    return {entry: outward_coefficient(entry) for entry in repaired_entries()}


def diagnostic() -> MaskedQuinticSliceRepair:
    affected = coefficient_one_dependent_entries()
    candidates = candidate_entries()
    repaired = repaired_entries()
    coefficients = tuple(outward_coefficient(entry) for entry in repaired)
    squared = tuple(sorted({float(squared_coefficient(entry)) for entry in repaired}))
    if len(candidates) != 54 or len(repaired) != 54:
        raise AssertionError(
            ("unexpected quintic slice inventory", len(candidates), len(repaired))
        )
    return MaskedQuinticSliceRepair(
        order=ORDER,
        dimension=ORDER**2,
        coefficient_one_dependent_entries=len(affected),
        candidate_entries=len(candidates),
        repaired_entries=len(repaired),
        remaining_quarantined_entries=len(affected) - len(repaired),
        repaired_profiles=len({profile for profile, _ in repaired}),
        bi_record_quintic_tail_bound=1 / (ORDER - 1) ** 2,
        distinct_squared_coefficients=squared,
        minimum_coefficient=min(coefficients),
        maximum_coefficient=max(coefficients),
        maximum_squared_coefficient=max(
            float(squared_coefficient(entry)) for entry in repaired
        ),
    )


def artifact_text(result: MaskedQuinticSliceRepair) -> str:
    payload = {
        "schema": "round4_q64_masked_quintic_slice_repair_v1",
        "result": asdict(result),
        "repaired_registry_entries": [
            {
                "profile": list(profile),
                "split": list(split),
                "squared_coefficient_numerator": squared_coefficient(
                    (profile, split)
                ).numerator,
                "squared_coefficient_denominator": squared_coefficient(
                    (profile, split)
                ).denominator,
                "outward_coefficient": outward_coefficient((profile, split)),
            }
            for profile, split in repaired_entries()
        ],
        "evidence_label": (
            "arbitrary-correlated-diagonal theorem for the actual masked kernel; "
            "exact rational quintic-singleton slices bound complete physical row "
            "or column energy, and displayed square-root coefficients are "
            "rounded upward"
        ),
        "remaining_open": (
            "306 coefficient-one-dependent entries remain quarantined; adaptive and "
            "interval Perron certification remain open"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-artifact", action="store_true")
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_masked_quintic_slice_repair.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 masked quintic slice repair: "
        f"candidates={result.candidate_entries},"
        f"repaired={result.repaired_entries},"
        f"remaining={result.remaining_quarantined_entries},"
        f"coefficient={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"maximum_squared={result.maximum_squared_coefficient:.12g},"
        "status=proved_actual_masked_arbitrary_diagonal_entries"
    )


if __name__ == "__main__":
    main()
