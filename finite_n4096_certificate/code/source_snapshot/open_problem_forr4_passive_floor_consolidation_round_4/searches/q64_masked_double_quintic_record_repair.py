#!/usr/bin/env python3
"""Repair all five masked double-quintic orbits by odd-record row energies."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from json import dumps
from math import comb, factorial, inf, nextafter, sqrt
from pathlib import Path

from q64_masked_double_quintic_endpoint_repair import candidate_entries


ROOT = Path(__file__).resolve().parents[1]
ORDER = 64
PROFILE = (1, 5, 5, 1)
TARGETS = (
    (PROFILE, (0, 2, 3, 1)),
    (PROFILE, (0, 2, 4, 0)),
    (PROFILE, (0, 3, 2, 1)),
    (PROFILE, (0, 3, 3, 0)),
    (PROFILE, (0, 4, 1, 1)),
)
RECORDS = (1, 3, 5)
ProfileSplit = tuple[tuple[int, ...], tuple[int, ...]]
Shape = tuple[tuple[int, int], ...]

# Exact maxima produced by ``enumerated_record_maxima``.  Keeping these as
# constants makes ordinary registry regeneration fast; the dedicated
# regression recomputes them from the finite shape enumeration.
Q64_VARIABLE_SINGLETON_FIXED_PAIR = {
    1: Fraction(4232, 63),
    3: Fraction(523388, 63),
    5: Fraction(10150400, 63),
}
Q64_FIXED_SINGLETON_FIXED_TRIPLE = {
    1: Fraction(151, 580608),
    3: Fraction(1811, 28672),
    5: Fraction(29905, 21504),
}
Q64_FIXED_SINGLETON_FIXED_ONE = {
    1: Fraction(5135, 4096),
    3: Fraction(114225, 2048),
    5: Fraction(247865, 256),
}
Q64_VARIABLE_SINGLETON_FIXED_FOUR_BOUNDS = {
    1: Fraction(2 * ORDER),
    3: Fraction(ORDER**2),
    5: Fraction(ORDER**2),
}
Q64_PAIR_MAXIMIZERS = {
    1: ((0, 0), (1, 0)),
    3: ((0, 0), (0, 1)),
    5: ((0, 0), (0, 1)),
}
Q64_TRIPLE_MAXIMIZERS = {
    1: ((0, 0), (1, 1), (1, 0)),
    3: ((0, 0), (0, 1), (0, 2)),
    5: ((0, 0), (0, 1), (0, 2)),
}


@dataclass(frozen=True)
class MaskedDoubleQuinticRecordRepair:
    order: int
    repaired_entries: int
    repaired_orbits: int
    remaining_double_quintic_entries: int
    central_records: tuple[int, ...]
    left_record_energies: tuple[str, ...]
    right_record_energies: tuple[str, ...]
    fixed_one_record_energies: tuple[str, ...]
    variable_four_record_bounds: tuple[str, ...]
    sector_squared_bounds: tuple[str, ...]
    distinct_squared_coefficient_bounds: tuple[str, ...]
    minimum_coefficient: float
    maximum_coefficient: float


def _set_partitions(size: int):
    """Yield restricted-growth strings for equality patterns."""

    values = [0] * size

    def visit(index: int, maximum: int):
        if index == size:
            yield tuple(values)
            return
        for value in range(maximum + 2):
            values[index] = value
            yield from visit(index + 1, max(maximum, value))

    if size == 0:
        yield ()
    else:
        yield from visit(1, 0)


def partial_shapes(selected: int) -> tuple[Shape, ...]:
    """Return all row/column equality types of a simple partial support.

    For at most three selected cells these are also the complete affine
    types over a binary vector space: any two distinct nonzero differences
    are linearly independent.
    """

    result = []
    for rows in _set_partitions(selected):
        for columns in _set_partitions(selected):
            shape = tuple(zip(rows, columns, strict=True))
            if len(set(shape)) == selected:
                result.append(shape)
    return tuple(result)


def fixed_one_record_energies(order: int) -> dict[int, Fraction]:
    """Return exact fixed-singleton, fixed-one quintic endpoint energies."""

    q = order
    if q < 4 or q & (q - 1):
        raise ValueError(("power-of-two order at least four required", q))
    w0 = Fraction(1, q**2)
    w1 = Fraction(1, q**2 * (q - 1) ** 2)
    w2 = Fraction(4, q**2 * (q - 1) ** 2 * (q - 2) ** 2)
    totals = {1: Fraction(0), 3: Fraction(0), 5: Fraction(0)}

    # Column pattern 5.
    totals[5] += q * comb(q, 5) * w0

    # Column pattern 4+1.
    column_choices = q * (q - 1)
    zero_xor_four_sets = q * (q - 1) * (q - 2) // 24
    four_sets = comb(q, 4)
    totals[3] += column_choices * (
        4 * zero_xor_four_sets * w0
        + 4 * (four_sets - zero_xor_four_sets) * w1
    )
    totals[5] += column_choices * (
        (q - 4) * zero_xor_four_sets * w0
        + (q - 4) * (four_sets - zero_xor_four_sets) * w1
    )

    # Column pattern 3+2, classified by row-set overlap.
    for overlap in range(3):
        count = (
            q
            * (q - 1)
            * comb(q, 3)
            * comb(3, overlap)
            * comb(q - 3, 2 - overlap)
        )
        totals[5 - 2 * overlap] += count * w1

    # Column pattern 2+2+1, classified by the two row-pair xors.
    column_choices = q * comb(q - 1, 2)
    pairs = comb(q, 2)
    total_counts = {
        1: pairs * q + pairs * 2 * (q - 2) * 2,
        3: pairs * 2 * (q - 2) ** 2 + pairs * comb(q - 2, 2) * 4,
        5: pairs * comb(q - 2, 2) * (q - 4),
    }
    equal_xor_counts = {
        1: pairs * q,
        3: pairs * (q // 2 - 1) * 4,
        5: pairs * (q // 2 - 1) * (q - 4),
    }
    for record in RECORDS:
        totals[record] += column_choices * (
            equal_xor_counts[record] * w1
            + (total_counts[record] - equal_xor_counts[record]) * w2
        )

    # Every quintic support contains five of the q^2 possible fixed cells.
    return {record: 5 * totals[record] / q**2 for record in RECORDS}


def variable_four_record_bounds(order: int) -> dict[int, Fraction]:
    """Bound a variable singleton endpoint with four fixed quintic cells.

    For central record one, a four-cell partial row pattern with two odd
    rows has at most ``2q`` completing cells.  If it has no odd rows, the
    pattern is ``4`` or ``2+2``.  Pattern ``4`` cannot also reach column
    record one because its four columns are distinct.  In pattern ``2+2``
    every completion leaves a size-two even row with nonzero xor, so its
    summed singleton weight is at most ``1/(q-1)^2`` and its total is below
    two.  Four odd rows cannot reach record one.  Thus ``B_1 <= 2q``.
    For records three and five, at most ``q^2`` cells complete the support
    and every summed singleton squared weight is at most one.
    """

    if order < 4 or order & (order - 1):
        raise ValueError(("power-of-two order at least four required", order))
    return {1: Fraction(2 * order), 3: Fraction(order**2), 5: Fraction(order**2)}


def parity_record(values: tuple[int, ...]) -> int:
    return sum(multiplicity % 2 for multiplicity in Counter(values).values())


def endpoint_square_weight(
    order: int,
    rows: tuple[int, ...],
    columns: tuple[int, ...],
    variable_singleton: bool,
) -> Fraction:
    """Squared singleton--quintic moment, summed over a singleton if free.

    One odd row is necessary.  Among the even row groups, let their nonzero
    column xors be ``x`` and possibly ``y``.  The structured permanent gives
    squared variable-singleton weights 1, 1/(q-1)^2, or
    4/((q-1)^2(q-2)^2); fixing the singleton divides by q^2.
    """

    groups: dict[int, list[int]] = defaultdict(list)
    for row, column in zip(rows, columns, strict=True):
        groups[row].append(column)
    if sum(len(values) % 2 for values in groups.values()) != 1:
        return Fraction(0)
    nonzero_xors = []
    for values in groups.values():
        if len(values) % 2 == 0:
            value = 0
            for column in values:
                value ^= column
            if value:
                nonzero_xors.append(value)
    if len(nonzero_xors) > 2:
        raise AssertionError((rows, columns, nonzero_xors))
    if not nonzero_xors:
        weight = Fraction(1)
    elif len(nonzero_xors) == 1 or nonzero_xors[0] == nonzero_xors[1]:
        weight = Fraction(1, (order - 1) ** 2)
    else:
        weight = Fraction(4, (order - 1) ** 2 * (order - 2) ** 2)
    return weight if variable_singleton else weight / order**2


def shape_record_energies(
    order: int,
    partial: Shape,
    variable_singleton: bool,
) -> dict[int, Fraction]:
    """Exactly sum all distinct quintic completions of one affine shape."""

    selected = len(partial)
    completion = 5 - selected
    fixed_rows = tuple(row for row, _ in partial)
    fixed_columns = tuple(column for _, column in partial)

    # The endpoint weight depends on rows only through their equality type.
    # Group the q^(5-k) actual row assignments before enumerating columns.
    row_patterns: dict[tuple[int, ...], int] = {}
    for extra_rows in product(range(order), repeat=completion):
        rows = fixed_rows + extra_rows
        labels: dict[int, int] = {}
        canonical = tuple(labels.setdefault(value, len(labels)) for value in rows)
        row_patterns[canonical] = row_patterns.get(canonical, 0) + 1

    totals: dict[int, Fraction] = defaultdict(Fraction)
    for rows, multiplicity in row_patterns.items():
        if parity_record(rows) != 1:
            continue
        subtotal: dict[int, Fraction] = defaultdict(Fraction)
        for extra_columns in product(range(order), repeat=completion):
            columns = fixed_columns + extra_columns
            if len(set(zip(rows, columns, strict=True))) != 5:
                continue
            record = parity_record(columns)
            subtotal[record] += endpoint_square_weight(
                order, rows, columns, variable_singleton
            )
        for record, value in subtotal.items():
            # Completions are unordered occurrence subsets.
            totals[record] += multiplicity * value / factorial(completion)
    return dict(totals)


@lru_cache(maxsize=None)
def enumerated_record_maxima(
    order: int,
    selected: int,
    variable_singleton: bool,
) -> tuple[tuple[int, Fraction, Shape], ...]:
    """Maximize exact record energies over every relevant affine shape."""

    if order < 4 or order & (order - 1):
        raise ValueError(("power-of-two order at least four required", order))
    rows = tuple(
        (shape, shape_record_energies(order, shape, variable_singleton))
        for shape in partial_shapes(selected)
    )
    return tuple(
        (
            record,
            *max(
                (energies.get(record, Fraction(0)), shape)
                for shape, energies in rows
            ),
        )
        for record in RECORDS
    )


def orbit(target: ProfileSplit) -> tuple[ProfileSplit, ...]:
    profile, split = target
    complement = tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )
    reverse_profile = tuple(reversed(profile))
    return tuple(
        sorted(
            {
                target,
                (profile, complement),
                (reverse_profile, tuple(reversed(split))),
                (reverse_profile, tuple(reversed(complement))),
            }
        )
    )


def repaired_entries() -> tuple[ProfileSplit, ...]:
    entries = tuple(sorted({entry for target in TARGETS for entry in orbit(target)}))
    if len(entries) != 12 or not set(entries).issubset(candidate_entries()):
        raise AssertionError(("double-quintic record orbit", entries))
    return entries


def endpoint_record_energies(
    selected: int,
    variable_singleton: bool,
    order: int = ORDER,
) -> dict[int, Fraction] | None:
    """Return an exact endpoint family or a stated rational upper bound."""

    if order == ORDER:
        pair = Q64_VARIABLE_SINGLETON_FIXED_PAIR
        triple = Q64_FIXED_SINGLETON_FIXED_TRIPLE
        fixed_one = Q64_FIXED_SINGLETON_FIXED_ONE
        variable_four = Q64_VARIABLE_SINGLETON_FIXED_FOUR_BOUNDS
    else:
        pair = {
            record: value
            for record, value, _ in enumerated_record_maxima(order, 2, True)
        }
        triple = {
            record: value
            for record, value, _ in enumerated_record_maxima(order, 3, False)
        }
        fixed_one = fixed_one_record_energies(order)
        variable_four = variable_four_record_bounds(order)
    if selected == 1 and not variable_singleton:
        return fixed_one
    if selected == 2:
        if variable_singleton:
            return pair
        return {record: value / order**2 for record, value in pair.items()}
    if selected == 3:
        if variable_singleton:
            return {
                record: value * order**2 for record, value in triple.items()
            }
        return triple
    if selected == 4 and variable_singleton:
        return variable_four
    return None


def row_sector_squared_bounds(
    entry: ProfileSplit,
    order: int = ORDER,
) -> dict[int, Fraction] | None:
    """Bound one complete physical row, or return None if not calibrated."""

    profile, split = entry
    if profile != PROFILE:
        raise ValueError(("not a double-quintic entry", entry))
    left = endpoint_record_energies(split[1], split[0] == 0, order)
    right = endpoint_record_energies(split[2], split[3] == 0, order)
    if left is None or right is None:
        return None
    return {
        record: left[record]
        * right[record]
        / comb(order, record) ** 2
        for record in RECORDS
    }


def row_squared_coefficient(
    entry: ProfileSplit,
    order: int = ORDER,
) -> Fraction | None:
    sectors = row_sector_squared_bounds(entry, order)
    return None if sectors is None else sum(sectors.values(), start=Fraction(0))


def squared_coefficient_bound(
    entry: ProfileSplit,
    order: int = ORDER,
) -> Fraction:
    """Take the better complete-row bound for the matrix or its transpose."""

    profile, split = entry
    complement = tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )
    candidates = tuple(
        value
        for value in (
            row_squared_coefficient(entry, order),
            row_squared_coefficient((profile, complement), order),
        )
        if value is not None
    )
    if not candidates:
        raise ValueError(("no calibrated row orientation", entry))
    return min(candidates)


def outward_coefficient(entry: ProfileSplit, order: int = ORDER) -> float:
    exact = squared_coefficient_bound(entry, order)
    value = sqrt(float(exact))
    while Fraction.from_float(value) ** 2 < exact:
        value = nextafter(value, inf)
    return value


def coefficient_map() -> dict[ProfileSplit, float]:
    return {entry: outward_coefficient(entry) for entry in repaired_entries()}


def diagnostic() -> MaskedDoubleQuinticRecordRepair:
    coefficients = coefficient_map()
    exact = tuple(sorted({squared_coefficient_bound(entry) for entry in repaired_entries()}))
    first_sectors = row_sector_squared_bounds(TARGETS[0])
    if first_sectors is None:
        raise AssertionError("missing first record sectors")
    return MaskedDoubleQuinticRecordRepair(
        order=ORDER,
        repaired_entries=len(repaired_entries()),
        repaired_orbits=len(TARGETS),
        remaining_double_quintic_entries=12 - len(repaired_entries()),
        central_records=RECORDS,
        left_record_energies=tuple(
            str(Q64_VARIABLE_SINGLETON_FIXED_PAIR[record]) for record in RECORDS
        ),
        right_record_energies=tuple(
            str(Q64_FIXED_SINGLETON_FIXED_TRIPLE[record]) for record in RECORDS
        ),
        fixed_one_record_energies=tuple(
            str(Q64_FIXED_SINGLETON_FIXED_ONE[record]) for record in RECORDS
        ),
        variable_four_record_bounds=tuple(
            str(Q64_VARIABLE_SINGLETON_FIXED_FOUR_BOUNDS[record])
            for record in RECORDS
        ),
        sector_squared_bounds=tuple(
            str(first_sectors[record]) for record in RECORDS
        ),
        distinct_squared_coefficient_bounds=tuple(map(str, exact)),
        minimum_coefficient=min(coefficients.values()),
        maximum_coefficient=max(coefficients.values()),
    )


def artifact_text(result: MaskedDoubleQuinticRecordRepair) -> str:
    payload = {
        "schema": "round4_q64_masked_double_quintic_record_repair_v2",
        "result": asdict(result),
        "repaired_registry_entries": [
            {
                "profile": list(profile),
                "split": list(split),
                "squared_coefficient_bound": str(
                    squared_coefficient_bound((profile, split))
                ),
                "outward_coefficient": outward_coefficient((profile, split)),
            }
            for profile, split in repaired_entries()
        ],
        "q64_maximizing_shapes": {
            "variable_singleton_fixed_pair": {
                str(record): [list(cell) for cell in Q64_PAIR_MAXIMIZERS[record]]
                for record in RECORDS
            },
            "fixed_singleton_fixed_triple": {
                str(record): [list(cell) for cell in Q64_TRIPLE_MAXIMIZERS[record]]
                for record in RECORDS
            },
        },
        "evidence_label": (
            "arbitrary-correlated-diagonal actual-mask theorem; exact finite "
            "affine-shape enumeration and parity bounds for complete endpoint row factors; "
            "odd-record sector decomposition; central quintic-quintic "
            "permanent bounded by 1/binomial(q,r); exact rational q64 row "
            "energies or rational upper bounds and outward-rounded coefficients"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--verify-exact-enumeration",
        action="store_true",
        help="recompute the q64 affine-shape maxima before writing",
    )
    arguments = parser.parse_args()
    if arguments.verify_exact_enumeration:
        pair = enumerated_record_maxima(ORDER, 2, True)
        triple = enumerated_record_maxima(ORDER, 3, False)
        if {r: v for r, v, _ in pair} != Q64_VARIABLE_SINGLETON_FIXED_PAIR:
            raise AssertionError(("q64 pair maxima", pair))
        if {r: v for r, v, _ in triple} != Q64_FIXED_SINGLETON_FIXED_TRIPLE:
            raise AssertionError(("q64 triple maxima", triple))
        if fixed_one_record_energies(ORDER) != Q64_FIXED_SINGLETON_FIXED_ONE:
            raise AssertionError("q64 fixed-one record energies")
        if variable_four_record_bounds(ORDER) != Q64_VARIABLE_SINGLETON_FIXED_FOUR_BOUNDS:
            raise AssertionError("q64 variable-four record bounds")
    result = diagnostic()
    if arguments.output is not None:
        arguments.output.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 masked double-quintic record repair: "
        f"repaired={result.repaired_entries},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"remaining_double_quintic={result.remaining_double_quintic_entries}"
    )


if __name__ == "__main__":
    main()
