#!/usr/bin/env python3
"""Actual-mask contractions for all eighty final q64 residual entries."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from fractions import Fraction
from json import dumps
from math import comb
from pathlib import Path

from q64_degree_ten_completion_row_insertion import orbit
from q64_masked_cubic_septimic_chain_repair import (
    degree_seven_endpoint_energy_bound,
    degree_seven_incidence,
)
from q64_masked_four_cubic_incidence_repair import (
    record_one_link_bound,
    record_three_link_bound,
)
from q64_masked_local_walsh_repair import (
    entry_integer_mask_factor,
    mechanism as local_walsh_mechanism,
)
from q64_masked_recovered_cubic_quintic_incidence_repair import (
    block_incidence_table,
    outward_sqrt,
)
from q64_same_side_whole_link_insertion import remaining_entries


ROOT = Path(__file__).resolve().parents[1]
ORDER = 64
ProfileSplit = tuple[tuple[int, ...], tuple[int, ...]]

TERMINAL_HIGH_TWO_CUBIC = ((1, 3, 3, 5), (0, 3, 0, 3))
DOUBLE_ENDPOINT_GENERATORS = (
    ((1, 3, 5, 1), (0, 3, 1, 1)),
    ((1, 3, 5, 1), (0, 3, 2, 0)),
    ((1, 3, 7, 1), (0, 3, 2, 1)),
    ((1, 3, 7, 1), (0, 3, 3, 0)),
)
ENDPOINT_CUBIC_HIGH = ((1, 3, 5, 3), (0, 3, 3, 0))
ENDPOINT_HIGH_CUBIC_COMPLETION = ((1, 5, 3, 3), (0, 3, 0, 3))
ENDPOINT_HIGH_CUBIC_ROW = ((1, 5, 3, 3), (0, 3, 3, 0))
CHAIN_GENERATORS = (
    TERMINAL_HIGH_TWO_CUBIC,
    *DOUBLE_ENDPOINT_GENERATORS,
    ENDPOINT_CUBIC_HIGH,
    ENDPOINT_HIGH_CUBIC_COMPLETION,
    ENDPOINT_HIGH_CUBIC_ROW,
)


@dataclass(frozen=True)
class FinalResidualChainContraction:
    order: int
    residual_entries: int
    residual_orbits: int
    local_walsh_entries: int
    local_walsh_orbits: int
    chain_entries: int
    chain_orbits: int
    terminal_high_two_cubic_orbits: int
    double_endpoint_orbits: int
    endpoint_cubic_high_orbits: int
    endpoint_high_cubic_orbits: int
    distinct_coefficients: tuple[float, ...]
    minimum_coefficient: float
    maximum_coefficient: float
    maximum_entry: ProfileSplit
    all_coefficients_strictly_below_one: bool


def validate_order(order: int) -> None:
    if order < 8 or order & (order - 1):
        raise ValueError(("power-of-two order at least eight required", order))


def canonical_orbit(entry: ProfileSplit) -> ProfileSplit:
    return min(orbit(entry))


def quintic_endpoint_fixed_three_energy(order: int = ORDER) -> Fraction:
    """Exact singleton--quintic squared slice through a fixed triple."""

    validate_order(order)
    q = order
    w0 = Fraction(1, q * q)
    w1 = Fraction(1, q * q * (q - 1) ** 2)
    return (
        comb(q - 3, 2) * w0
        + q * (q - 1) * (w0 + (q - 4) * w1)
        + (q - 1) * comb(q, 2) * w1
    )


def septimic_record_three_zero_active_total(order: int = ORDER) -> int:
    """Count all record-(3,1) zero-active septimic supports."""

    validate_order(order)
    q = order
    zero_xor_four_sets = q * (q - 1) * (q - 2) // 24
    return q * (q - 1) * zero_xor_four_sets * 6 * (q - 4)


def septimic_record_three_zero_active_fixed_cell(
    order: int = ORDER,
) -> int:
    """Exact fixed-cell incidence of the translation-invariant zero family."""

    total = septimic_record_three_zero_active_total(order)
    numerator = 7 * total
    denominator = order * order
    if numerator % denominator:
        raise AssertionError("nonintegral fixed-cell zero-active incidence")
    return numerator // denominator


def septimic_record_three_endpoint_energy_bound(
    selected: int,
    order: int = ORDER,
) -> Fraction:
    """Endpoint energy for a record-(3,1) septimic fixed triple.

    With no active even endpoint group the only physical pattern is 4+3:
    the four-neighborhood has xor zero and meets the three-neighborhood in
    two labels. Every fixed nonempty partial is contained in at most the
    fixed-cell incidence of this family. Every other support has a nonzero
    injective-Walsh residual bounded by 1/(q-1).
    """

    validate_order(order)
    if selected != 3:
        raise ValueError(("only the fixed-three bound is used", selected))
    return Fraction(septimic_record_three_zero_active_fixed_cell(order)) + Fraction(
        degree_seven_incidence(selected, 3, order),
        (order - 1) ** 2,
    )


def terminal_high_two_cubic_squared(order: int = ORDER) -> Fraction:
    """Complete row for singleton--cubic--cubic--split-quintic."""

    validate_order(order)
    q = order
    b = record_one_link_bound(q)
    g = record_three_link_bound(q)
    endpoint_one = Fraction(1, q - 1)
    fixed = 3

    def cubic_count(left_record: int, right_record: int) -> int:
        return block_incidence_table(
            3, left_record, right_record, q
        )[0]

    def quintic_count(record: int) -> int:
        return block_incidence_table(5, record, None, q)[fixed]

    result = (
        cubic_count(1, 1)
        * quintic_count(1)
        * (endpoint_one * b * b) ** 2
    )
    result += (
        cubic_count(1, 3)
        * quintic_count(3)
        * (endpoint_one * b * g) ** 2
    )
    result += (
        cubic_count(3, 3)
        * quintic_count(3)
        * g**4
    )

    record_three_one = cubic_count(3, 1)
    no_even_cubic = q * comb(q, 3)
    if no_even_cubic > record_three_one:
        raise AssertionError("invalid cubic shape split")
    result += g**2 * (
        no_even_cubic * quintic_endpoint_fixed_three_energy(q)
        + (record_three_one - no_even_cubic)
        * quintic_count(1)
        * b**2
    )
    return result


def double_endpoint_squared(
    degree: int,
    selected: int,
    terminal_singleton_selected: int,
    order: int = ORDER,
) -> Fraction:
    """Complete row for singleton--cubic--split-higher--singleton."""

    validate_order(order)
    q = order
    if degree not in (5, 7):
        raise ValueError(("quintic or septimic required", degree))
    if terminal_singleton_selected not in (0, 1):
        raise ValueError(("singleton split", terminal_singleton_selected))
    result = Fraction(0)
    for record in (1, 3):
        endpoint_amplitude = Fraction(1, q - 1) if record == 1 else Fraction(1)
        if record == 1:
            middle = record_one_link_bound(q) if degree == 5 else Fraction(1, q)
        else:
            middle = record_three_link_bound(q)
        if terminal_singleton_selected == 1:
            terminal_energy = Fraction(
                block_incidence_table(degree, record, 1, q)[selected],
                q * q,
            )
        elif degree == 5:
            terminal_energy = Fraction(
                block_incidence_table(degree, record, 1, q)[selected]
            )
        elif record == 1:
            terminal_energy = degree_seven_endpoint_energy_bound(
                selected, record, q
            )
        else:
            terminal_energy = septimic_record_three_endpoint_energy_bound(
                selected, q
            )
        result += endpoint_amplitude**2 * middle**2 * terminal_energy
    return result


def endpoint_cubic_high_squared(order: int = ORDER) -> Fraction:
    """Masked endpoint--cubic--quintic row with the last link completed."""

    validate_order(order)
    q = order
    b = record_one_link_bound(q)
    g = record_three_link_bound(q)
    return (
        Fraction(block_incidence_table(5, 1, None, q)[3])
        * (b / (q - 1)) ** 2
        + Fraction(block_incidence_table(5, 3, None, q)[3]) * g**2
    )


def endpoint_high_cubic_completion_squared(order: int = ORDER) -> Fraction:
    """Fixed endpoint and cubic with the complementary quintic pair fixed."""

    validate_order(order)
    q = order
    return sum(
        Fraction(
            block_incidence_table(5, 1, record, q)[2],
            q * q * comb(q, record) ** 2,
        )
        for record in (1, 3)
    )


def endpoint_high_cubic_row_squared(order: int = ORDER) -> Fraction:
    """Endpoint--quintic--cubic complete row through a fixed triple."""

    validate_order(order)
    return quintic_endpoint_fixed_three_energy(order) / (order * order)


def chain_mechanism(entry: ProfileSplit) -> str:
    canonical = canonical_orbit(entry)
    if canonical == TERMINAL_HIGH_TWO_CUBIC:
        return "terminal_high_two_cubic_row"
    if canonical in DOUBLE_ENDPOINT_GENERATORS:
        return "double_endpoint_high_row"
    if canonical == ENDPOINT_CUBIC_HIGH:
        return "endpoint_cubic_high_row"
    if canonical == ENDPOINT_HIGH_CUBIC_COMPLETION:
        return "endpoint_high_cubic_completion"
    if canonical == ENDPOINT_HIGH_CUBIC_ROW:
        return "endpoint_high_cubic_row"
    raise ValueError(("not a residual chain orbit", entry))


def chain_squared_coefficient(
    entry: ProfileSplit,
    order: int = ORDER,
) -> Fraction:
    canonical = canonical_orbit(entry)
    if canonical == TERMINAL_HIGH_TWO_CUBIC:
        return terminal_high_two_cubic_squared(order)
    if canonical in DOUBLE_ENDPOINT_GENERATORS:
        profile, split = canonical
        return double_endpoint_squared(
            profile[2], split[2], split[3], order
        )
    if canonical == ENDPOINT_CUBIC_HIGH:
        return endpoint_cubic_high_squared(order)
    if canonical == ENDPOINT_HIGH_CUBIC_COMPLETION:
        return endpoint_high_cubic_completion_squared(order)
    if canonical == ENDPOINT_HIGH_CUBIC_ROW:
        return endpoint_high_cubic_row_squared(order)
    raise ValueError(("not a residual chain orbit", entry))


def mechanism(entry: ProfileSplit) -> str:
    local = local_walsh_mechanism(entry)
    if local is not None:
        return f"local_walsh_{local}"
    return chain_mechanism(entry)


def squared_coefficient(
    entry: ProfileSplit,
    order: int = ORDER,
) -> Fraction:
    if entry not in set(remaining_entries()):
        raise ValueError(("not a final residual entry", entry))
    if local_walsh_mechanism(entry) is not None:
        value = Fraction(entry_integer_mask_factor(entry), order)
        return value**2
    return chain_squared_coefficient(entry, order)


def coefficient(entry: ProfileSplit, order: int = ORDER) -> float:
    return outward_sqrt(squared_coefficient(entry, order))


def repaired_entries(order: int = ORDER) -> tuple[ProfileSplit, ...]:
    entries = tuple(sorted(remaining_entries()))
    if len(entries) != 80:
        raise AssertionError(("final residual inventory", len(entries)))
    return tuple(entry for entry in entries if coefficient(entry, order) <= 1)


def coefficient_map(order: int = ORDER) -> dict[ProfileSplit, float]:
    return {entry: coefficient(entry, order) for entry in repaired_entries(order)}


def diagnostic() -> FinalResidualChainContraction:
    entries = tuple(sorted(remaining_entries()))
    repaired = repaired_entries()
    if len(repaired) != len(entries):
        raise AssertionError(("unrepaired final residual", len(repaired), len(entries)))
    coefficients = coefficient_map()
    local = tuple(entry for entry in entries if local_walsh_mechanism(entry) is not None)
    chain = tuple(entry for entry in entries if local_walsh_mechanism(entry) is None)
    mechanisms = Counter(chain_mechanism(entry) for entry in chain)
    maximum_entry = max(coefficients, key=coefficients.get)
    values = tuple(sorted(set(coefficients.values())))
    return FinalResidualChainContraction(
        order=ORDER,
        residual_entries=len(entries),
        residual_orbits=len({canonical_orbit(entry) for entry in entries}),
        local_walsh_entries=len(local),
        local_walsh_orbits=len({canonical_orbit(entry) for entry in local}),
        chain_entries=len(chain),
        chain_orbits=len({canonical_orbit(entry) for entry in chain}),
        terminal_high_two_cubic_orbits=(
            mechanisms["terminal_high_two_cubic_row"] // 4
        ),
        double_endpoint_orbits=mechanisms["double_endpoint_high_row"] // 4,
        endpoint_cubic_high_orbits=mechanisms["endpoint_cubic_high_row"] // 4,
        endpoint_high_cubic_orbits=(
            mechanisms["endpoint_high_cubic_completion"]
            + mechanisms["endpoint_high_cubic_row"]
        )
        // 4,
        distinct_coefficients=values,
        minimum_coefficient=min(values),
        maximum_coefficient=coefficients[maximum_entry],
        maximum_entry=maximum_entry,
        all_coefficients_strictly_below_one=max(values) < 1,
    )


def artifact_text(result: FinalResidualChainContraction) -> str:
    payload = {
        "schema": "round4_q64_final_residual_chain_contraction_v1",
        "result": asdict(result),
        "registry_entries": [
            {
                "profile": list(profile),
                "split": list(split),
                "canonical_profile": list(canonical_orbit((profile, split))[0]),
                "canonical_split": list(canonical_orbit((profile, split))[1]),
                "mechanism": mechanism((profile, split)),
                "squared_coefficient_numerator": squared_coefficient(
                    (profile, split)
                ).numerator,
                "squared_coefficient_denominator": squared_coefficient(
                    (profile, split)
                ).denominator,
                "outward_coefficient": coefficient((profile, split)),
            }
            for profile, split in repaired_entries()
        ],
        "evidence_label": (
            "arbitrary-correlated-diagonal one-batch theorem for all eighty "
            "final residual entries; forty-eight entries use the exact local "
            "Walsh q^-1 mechanism with their sole physical mask restored, "
            "and thirty-two entries use exact complete-chain row energies "
            "with physical incidence tables and active-xor endpoint bounds; "
            "all squared coefficients are rational and only final square "
            "roots are rounded outward"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = diagnostic()
    text = artifact_text(result)
    if arguments.output is not None:
        arguments.output.write_text(text, encoding="utf-8")
    print(
        "q64 final residual chain contraction: "
        f"entries={result.residual_entries},orbits={result.residual_orbits},"
        f"local={result.local_walsh_entries}/{result.local_walsh_orbits},"
        f"chain={result.chain_entries}/{result.chain_orbits},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"maximum_entry={result.maximum_entry},"
        "status=all_final_residual_entries_proved"
    )


if __name__ == "__main__":
    main()
