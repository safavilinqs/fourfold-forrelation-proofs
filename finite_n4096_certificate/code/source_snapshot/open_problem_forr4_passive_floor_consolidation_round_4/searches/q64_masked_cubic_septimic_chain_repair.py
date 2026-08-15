#!/usr/bin/env python3
"""Repair all residual cubic--septimic entries by masked chain rows."""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict, dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from json import dumps
from math import comb, factorial, inf, nextafter, prod, sqrt
from pathlib import Path

from q64_masked_four_cubic_incidence_repair import (
    coefficient_one_dependent_entries,
    previous_repaired_entries as pre_four_cubic_repaired_entries,
    repaired_entries as four_cubic_repaired_entries,
)


ROOT = Path(__file__).resolve().parents[1]
ORDER = 64
CANONICAL_PROFILE = (1, 3, 7, 1)
RECORDS = (1, 3)
ProfileSplit = tuple[tuple[int, ...], tuple[int, ...]]
Partition = tuple[int, ...]


@dataclass(frozen=True)
class MaskedCubicSeptimicChainRepair:
    order: int
    candidate_entries: int
    repaired_entries: int
    repaired_orbits: int
    degree_seven_record_one_one_fixed_three_incidence: int
    degree_seven_record_one_one_fixed_four_incidence: int
    degree_seven_record_three_one_fixed_three_incidence: int
    degree_seven_record_three_one_fixed_four_incidence: int
    record_one_zero_active_completion_bound: int
    record_one_fixed_three_endpoint_energy_bound: float
    record_one_fixed_four_endpoint_energy_bound: float
    cubic_fixed_one_record_one_energy: float
    cubic_fixed_one_record_three_energy: float
    minimum_coefficient: float
    maximum_coefficient: float
    remaining_quarantined_entries: int


def canonical(values: tuple[int, ...]) -> tuple[int, ...]:
    labels: dict[int, int] = {}
    return tuple(labels.setdefault(value, len(labels)) for value in values)


def set_partitions(size: int):
    """Yield restricted-growth encodings of all set partitions."""

    if size == 0:
        yield ()
        return
    values = [0] * size

    def visit(index: int, maximum: int):
        if index == size:
            yield tuple(values)
            return
        for value in range(maximum + 2):
            values[index] = value
            yield from visit(index + 1, max(maximum, value))

    values[0] = 0
    yield from visit(1, 0)


def parity_record(partition: Partition) -> int:
    counts: dict[int, int] = defaultdict(int)
    for value in partition:
        counts[value] += 1
    return sum(count % 2 for count in counts.values())


def simple_bipartite_shape(rows: Partition, columns: Partition) -> bool:
    return len(set(zip(rows, columns, strict=True))) == len(rows)


def falling(total: int, selected: int) -> int:
    return prod(range(total - selected + 1, total + 1))


@lru_cache(maxsize=None)
def degree_seven_bidegree_incidence_table(
    order: int,
    row_record: int,
    column_record: int,
) -> tuple[int, ...]:
    """Exact maximum physical completion incidences for degree seven.

    A full support is a simple bipartite graph with seven labeled edges.
    Row and column set partitions encode its equality shape.  Restricting the
    first ``k`` edges gives the fixed occurrence partial support.  Falling
    factorials count embeddings of the new vertices, and division by
    ``(7-k)!`` forgets the completion-edge ordering.
    """

    degree = 7
    grouped: dict[int, list[Partition]] = defaultdict(list)
    for partition in set_partitions(degree):
        grouped[parity_record(partition)].append(partition)
    accumulators = [defaultdict(int) for _ in range(degree + 1)]
    for rows in grouped[row_record]:
        total_rows = len(set(rows))
        for columns in grouped[column_record]:
            if not simple_bipartite_shape(rows, columns):
                continue
            total_columns = len(set(columns))
            for selected in range(degree + 1):
                partial_rows = canonical(rows[:selected])
                partial_columns = canonical(columns[:selected])
                old_rows = len(set(rows[:selected]))
                old_columns = len(set(columns[:selected]))
                accumulators[selected][(partial_rows, partial_columns)] += (
                    falling(order - old_rows, total_rows - old_rows)
                    * falling(
                        order - old_columns,
                        total_columns - old_columns,
                    )
                )
    result = []
    for selected, values in enumerate(accumulators):
        divisor = factorial(degree - selected)
        if any(value % divisor for value in values.values()):
            raise AssertionError("nonintegral degree-seven incidence")
        result.append(max(value // divisor for value in values.values()))
    return tuple(result)


def degree_seven_incidence(
    selected: int,
    central_record: int,
    order: int = ORDER,
) -> int:
    """Completion incidence with central record r and endpoint record one."""

    return degree_seven_bidegree_incidence_table(
        order, central_record, 1
    )[selected]


def zero_active_completion_bound(order: int = ORDER) -> int:
    """Bound record-one endpoint completions with no active even group.

    A seven-cell support with endpoint-axis record one and no active even
    endpoint group can only have endpoint degrees ``3+4``.  The four-cell
    neighborhood has xor zero, and central record one forces the three-cell
    neighborhood to be a subset of it.  A fixed three- or four-cell partial
    support is contained in at most ``5(q-1)`` such supports: if it uses one
    endpoint coordinate, the two possible degree roles contribute at most
    ``(q-1)+4(q-1)``; if it uses both, the zero-xor four-set count is smaller.
    """

    if order < 4 or order & (order - 1):
        raise ValueError(("power-of-two order at least four required", order))
    return 5 * (order - 1)


def nonzero_active_character_square_bound(order: int = ORDER) -> Fraction:
    """Uniform squared endpoint residual once an even group is active."""

    q = order
    if q < 4 or q & (q - 1):
        raise ValueError(("power-of-two order at least four required", q))
    return max(
        Fraction(1, (q - 1) ** 2),
        Fraction(9, (q - 1) ** 2 * (q - 3) ** 2),
    )


def degree_seven_endpoint_energy_bound(
    selected: int,
    central_record: int,
    order: int = ORDER,
) -> Fraction:
    """Summed variable-singleton endpoint energy, mask retained.

    In central record one, separate the at most ``5(q-1)`` zero-active
    completions.  Every other completion gains the squared injective-Walsh
    residual returned above.  Record three keeps the safe incidence bound;
    its central moment is already suppressed by ``binom(q,3)``.
    """

    incidence = degree_seven_incidence(selected, central_record, order)
    if central_record == 1:
        return Fraction(zero_active_completion_bound(order)) + (
            nonzero_active_character_square_bound(order) * incidence
        )
    if central_record == 3:
        return Fraction(incidence)
    raise ValueError(("unsupported central record", central_record))


def cubic_endpoint_slice_energies(order: int = ORDER) -> tuple[Fraction, ...]:
    q = order
    return (
        Fraction(q * q + 2, 6),
        Fraction(q * q + 2, 2 * q * q),
        Fraction(q * q - 2 * q + 2, q * q * (q - 1)),
        Fraction(1, q * q),
    )


def cubic_fixed_one_record_energies(
    order: int = ORDER,
) -> dict[int, Fraction]:
    """Exact fixed-one cubic endpoint energy by the following-link record."""

    q = order
    return {
        1: Fraction(3, q * q),
        3: Fraction(q * q - 4, 2 * q * q),
    }


def outward_sqrt(exact: Fraction) -> float:
    value = sqrt(float(exact))
    while Fraction.from_float(value) ** 2 < exact:
        value = nextafter(value, inf)
    return value


def outward_sum(values) -> float:
    result = 0.0
    for value in values:
        result = nextafter(result + value, inf)
    return result


def sector_squared_coefficient(
    endpoint_energy: Fraction,
    septimic_selected: int,
    central_record: int,
    endpoint_singleton_selected: int,
    order: int = ORDER,
) -> Fraction:
    q = order
    if endpoint_singleton_selected == 0:
        endpoint_energy_bound = degree_seven_endpoint_energy_bound(
            septimic_selected, central_record, q
        )
    elif endpoint_singleton_selected == 1:
        endpoint_energy_bound = Fraction(
            degree_seven_incidence(septimic_selected, central_record, q),
            q * q,
        )
    else:
        raise ValueError(("invalid singleton split", endpoint_singleton_selected))
    return (
        endpoint_energy
        * endpoint_energy_bound
        / comb(q, central_record) ** 2
    )


def orientation_coefficient(
    split: tuple[int, ...],
    order: int = ORDER,
) -> tuple[float, str]:
    """Bound one canonical row orientation, masks included."""

    (
        singleton_selected,
        cubic_selected,
        septimic_selected,
        endpoint_singleton_selected,
    ) = split
    slices = cubic_endpoint_slice_energies(order)
    total_endpoint_energy = (
        order ** (2 * (1 - singleton_selected))
        * slices[cubic_selected]
    )
    coarse = outward_sum(
        outward_sqrt(
            sector_squared_coefficient(
                total_endpoint_energy,
                septimic_selected,
                record,
                endpoint_singleton_selected,
                order,
            )
        )
        for record in RECORDS
    )
    candidates = [(coarse, "complete_cubic_endpoint")]
    if singleton_selected == 1 and cubic_selected == 1:
        energies = cubic_fixed_one_record_energies(order)
        refined = outward_sum(
            outward_sqrt(
                sector_squared_coefficient(
                    energies[record],
                    septimic_selected,
                    record,
                    endpoint_singleton_selected,
                    order,
                )
            )
            for record in RECORDS
        )
        candidates.append((refined, "fixed_one_record_split"))
    return min(candidates)


def canonical_split(entry: ProfileSplit) -> tuple[int, ...]:
    profile, split = entry
    if profile == CANONICAL_PROFILE:
        return split
    if tuple(reversed(profile)) == CANONICAL_PROFILE:
        return tuple(reversed(split))
    raise ValueError(("not a cubic-septimic endpoint profile", entry))


def coefficient_with_mechanism(
    entry: ProfileSplit,
    order: int = ORDER,
) -> tuple[float, str, tuple[int, ...]]:
    split = canonical_split(entry)
    complement = tuple(
        degree - selected
        for degree, selected in zip(CANONICAL_PROFILE, split, strict=True)
    )
    candidates = []
    for label, orientation in (("row", split), ("column", complement)):
        value, mechanism = orientation_coefficient(orientation, order)
        candidates.append((value, f"{label}_{mechanism}", orientation))
    return min(candidates)


def coefficient(entry: ProfileSplit, order: int = ORDER) -> float:
    return coefficient_with_mechanism(entry, order)[0]


def previous_repaired_entries() -> frozenset[ProfileSplit]:
    return frozenset(pre_four_cubic_repaired_entries()) | frozenset(
        four_cubic_repaired_entries()
    )


def candidate_entries() -> tuple[ProfileSplit, ...]:
    previous = previous_repaired_entries()
    return tuple(
        sorted(
            entry
            for entry in coefficient_one_dependent_entries()
            if entry not in previous
            and sorted(entry[0]) == [1, 1, 3, 7]
        )
    )


def repaired_entries() -> tuple[ProfileSplit, ...]:
    return tuple(entry for entry in candidate_entries() if coefficient(entry) <= 1)


def coefficient_map() -> dict[ProfileSplit, float]:
    return {entry: coefficient(entry) for entry in repaired_entries()}


def diagnostic() -> MaskedCubicSeptimicChainRepair:
    candidates = candidate_entries()
    repaired = repaired_entries()
    if len(candidates) != 12 or len(repaired) != 12:
        raise AssertionError(("cubic-septimic inventory", len(candidates), len(repaired)))
    values = tuple(coefficient(entry) for entry in repaired)
    fixed_one = cubic_fixed_one_record_energies()
    return MaskedCubicSeptimicChainRepair(
        order=ORDER,
        candidate_entries=len(candidates),
        repaired_entries=len(repaired),
        repaired_orbits=3,
        degree_seven_record_one_one_fixed_three_incidence=(
            degree_seven_incidence(3, 1)
        ),
        degree_seven_record_one_one_fixed_four_incidence=(
            degree_seven_incidence(4, 1)
        ),
        degree_seven_record_three_one_fixed_three_incidence=(
            degree_seven_incidence(3, 3)
        ),
        degree_seven_record_three_one_fixed_four_incidence=(
            degree_seven_incidence(4, 3)
        ),
        record_one_zero_active_completion_bound=(
            zero_active_completion_bound()
        ),
        record_one_fixed_three_endpoint_energy_bound=float(
            degree_seven_endpoint_energy_bound(3, 1)
        ),
        record_one_fixed_four_endpoint_energy_bound=float(
            degree_seven_endpoint_energy_bound(4, 1)
        ),
        cubic_fixed_one_record_one_energy=float(fixed_one[1]),
        cubic_fixed_one_record_three_energy=float(fixed_one[3]),
        minimum_coefficient=min(values),
        maximum_coefficient=max(values),
        remaining_quarantined_entries=(
            len(coefficient_one_dependent_entries())
            - len(previous_repaired_entries())
            - len(repaired)
        ),
    )


def artifact_text(result: MaskedCubicSeptimicChainRepair) -> str:
    payload = {
        "schema": "round4_q64_masked_cubic_septimic_chain_repair_v2",
        "result": asdict(result),
        "incidence_tables": {
            f"record_{record}_one": list(
                degree_seven_bidegree_incidence_table(ORDER, record, 1)
            )
            for record in RECORDS
        },
        "repaired_registry_entries": [
            {
                "profile": list(profile),
                "split": list(split),
                "outward_coefficient": coefficient((profile, split)),
                "mechanism": coefficient_with_mechanism((profile, split))[1],
                "selected_orientation": list(
                    coefficient_with_mechanism((profile, split))[2]
                ),
            }
            for profile, split in repaired_entries()
        ],
        "evidence_label": (
            "arbitrary-correlated-diagonal theorem for all twelve residual "
            "cubic-septimic entries; exact physical cubic endpoint slices, "
            "a record-compatible fixed-one refinement, exact degree-seven "
            "bidegree incidences, and a zero-active/injective-Walsh endpoint "
            "split retain both occurrence masks"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.output is not None:
        arguments.output.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 masked cubic-septimic chain repair: "
        f"repaired={result.repaired_entries},"
        f"orbits={result.repaired_orbits},"
        f"coefficient_range={result.minimum_coefficient:.12g}/"
        f"{result.maximum_coefficient:.12g},"
        f"remaining={result.remaining_quarantined_entries}"
    )


if __name__ == "__main__":
    main()
