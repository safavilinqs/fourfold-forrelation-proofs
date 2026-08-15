#!/usr/bin/env python3
"""Independent physical audit of the recovered cubic--quintic repair."""

from __future__ import annotations

from collections import Counter, defaultdict
from decimal import Decimal, getcontext
from fractions import Fraction
from itertools import combinations, permutations, product
from json import loads
from math import comb, factorial
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from q64_masked_recovered_cubic_quintic_incidence_repair import (  # noqa: E402
    block_incidence_table,
    candidate_entries,
    canonical_entry,
    endpoint_energy,
    endpoint_remaining_link_bound,
    record_sectors,
    repaired_entries,
    rejected_chain_entries,
)
from q64_shared_quintic_row_chain_insertion import (  # noqa: E402
    odd_record_incidence,
    two_axis_relaxed_incidence,
)


Q4 = 4
Q8 = 8
ORBIT_REPRESENTATIVES = (
    ((1, 3, 3, 5), (0, 1, 2, 3)),
    ((1, 3, 3, 5), (0, 2, 1, 3)),
    ((1, 3, 3, 5), (0, 2, 2, 2)),
    ((1, 3, 5, 3), (0, 1, 3, 2)),
    ((1, 3, 5, 3), (0, 2, 2, 2)),
    ((1, 3, 5, 3), (0, 2, 3, 1)),
    ((1, 5, 3, 3), (0, 2, 2, 2)),
    ((1, 5, 3, 3), (0, 3, 1, 2)),
    ((1, 5, 3, 3), (0, 3, 2, 1)),
    ((1, 5, 3, 3), (0, 4, 1, 1)),
)


def walsh(row: int, column: int) -> int:
    return -1 if (row & column).bit_count() % 2 else 1


def parity_record(support: tuple[int, ...], order: int, axis: int) -> int:
    counts = Counter(divmod(coordinate, order)[axis] for coordinate in support)
    return sum(value % 2 for value in counts.values())


def support_records(support: tuple[int, ...], order: int) -> tuple[int, int]:
    return parity_record(support, order, 0), parity_record(support, order, 1)


def transpose_support(support: tuple[int, ...], order: int) -> tuple[int, ...]:
    return tuple(
        sorted(
            column * order + row
            for row, column in map(lambda value: divmod(value, order), support)
        )
    )


def direct_permutation_moment(
    order: int,
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> Fraction:
    """Average the physical link directly over permutations and signs.

    The sign average is performed exactly by rejecting a permutation unless
    every hidden sign occurs an even number of times.  No permanent formula
    or production link-moment helper is used.
    """

    total = 0
    for permutation in permutations(range(order)):
        inverse = [0] * order
        for domain, codomain in enumerate(permutation):
            inverse[codomain] = domain
        sign_parity = [0] * order
        value = 1
        for coordinate in left:
            row, column = divmod(coordinate, order)
            hidden = inverse[column]
            sign_parity[hidden] ^= 1
            value *= walsh(row, hidden)
        for coordinate in right:
            row, column = divmod(coordinate, order)
            sign_parity[row] ^= 1
            value *= walsh(permutation[row], column)
        if not any(sign_parity):
            total += value
    return Fraction(total, factorial(order))


class DirectQ4Plant:
    """Small exact physical plant built by full signed-permutation averaging."""

    def __init__(self) -> None:
        q = Q4
        dimension = q * q
        hadamard = np.array(
            [[walsh(row, column) for column in range(q)] for row in range(q)],
            dtype=np.int16,
        )
        left_values = []
        right_values = []
        for permutation in permutations(range(q)):
            for signs in product((-1, 1), repeat=q):
                signed_permutation = np.zeros((q, q), dtype=np.int16)
                for row, column in enumerate(permutation):
                    signed_permutation[row, column] = signs[row]
                left_values.append((hadamard @ signed_permutation).reshape(dimension))
                right_values.append((signed_permutation @ hadamard).reshape(dimension))
        self.left_values = np.asarray(left_values, dtype=np.int16)
        self.right_values = np.asarray(right_values, dtype=np.int16)
        self.group_size = len(left_values)
        self.supports = {
            degree: tuple(combinations(range(dimension), degree))
            for degree in (1, 3, 5)
        }
        self.support_index = {
            degree: {support: index for index, support in enumerate(supports)}
            for degree, supports in self.supports.items()
        }
        self._features: dict[tuple[str, int], np.ndarray] = {}
        self._moments: dict[tuple[int, int], np.ndarray] = {}

    def features(self, side: str, degree: int) -> np.ndarray:
        key = side, degree
        if key not in self._features:
            values = self.left_values if side == "left" else self.right_values
            indices = np.asarray(self.supports[degree], dtype=np.int16)
            # At q=4 the largest feature array has 1,677,312 int16 entries,
            # below 3.4 MB.  The largest moment matrix is 560 by 4,368 int16,
            # below 5 MB.
            self._features[key] = np.prod(
                values[:, indices], axis=2, dtype=np.int16
            ).T
        return self._features[key]

    def moments(self, left_degree: int, right_degree: int) -> np.ndarray:
        key = left_degree, right_degree
        if key not in self._moments:
            self._moments[key] = (
                self.features("left", left_degree)
                @ self.features("right", right_degree).T
            )
        return self._moments[key]

    def moment(
        self,
        left: tuple[int, ...],
        right: tuple[int, ...],
    ) -> Fraction:
        matrix = self.moments(len(left), len(right))
        numerator = int(
            matrix[
                self.support_index[len(left)][left],
                self.support_index[len(right)][right],
            ]
        )
        return Fraction(numerator, self.group_size)

    def filtered_maximum(
        self,
        left_degree: int,
        right_degree: int,
        left_records: tuple[int | None, int | None],
        right_records: tuple[int | None, int | None],
    ) -> Fraction:
        def admissible(
            support: tuple[int, ...], records: tuple[int | None, int | None]
        ) -> bool:
            observed = support_records(support, Q4)
            return all(
                expected is None or expected == value
                for expected, value in zip(records, observed, strict=True)
            )

        left_indices = [
            index
            for index, support in enumerate(self.supports[left_degree])
            if admissible(support, left_records)
        ]
        right_indices = [
            index
            for index, support in enumerate(self.supports[right_degree])
            if admissible(support, right_records)
        ]
        values = np.abs(
            self.moments(left_degree, right_degree)[
                np.ix_(left_indices, right_indices)
            ]
        )
        return Fraction(int(values.max(initial=0)), self.group_size)


def direct_endpoint_energy(
    plant: DirectQ4Plant,
    degree: int,
    selected: int,
    singleton_selected: int,
    outgoing_record: int,
) -> Fraction:
    """Maximum complete physical endpoint row energy at q=4."""

    accumulators: defaultdict[tuple[tuple[int, ...], tuple[int, ...]], int]
    accumulators = defaultdict(int)
    moments = plant.moments(1, degree)
    for higher_index, higher in enumerate(plant.supports[degree]):
        if support_records(higher, Q4) != (1, outgoing_record):
            continue
        for singleton_index, singleton in enumerate(plant.supports[1]):
            numerator = int(moments[singleton_index, higher_index])
            if not numerator:
                continue
            for partial in combinations(higher, selected):
                row_key = (
                    singleton if singleton_selected else (),
                    partial,
                )
                accumulators[row_key] += numerator * numerator
    return Fraction(max(accumulators.values(), default=0), plant.group_size**2)


def direct_incidence_maxima(
    order: int,
    degree: int,
    left_record: int | None,
    right_record: int | None,
) -> tuple[int, ...]:
    """Enumerate simple supports and their fixed-subset completion rows."""

    supports = tuple(combinations(range(order * order), degree))
    result = []
    for selected in range(degree + 1):
        counts: Counter[tuple[int, ...]] = Counter()
        for support in supports:
            row_record, column_record = support_records(support, order)
            if left_record is not None and row_record != left_record:
                continue
            if right_record is not None and column_record != right_record:
                continue
            counts.update(combinations(support, selected))
        result.append(max(counts.values(), default=0))
    return tuple(result)


def direct_one_axis_incidence_maxima(
    order: int, degree: int, record: int
) -> tuple[int, ...]:
    return direct_incidence_maxima(order, degree, record, None)


def orbit(entry: tuple[tuple[int, ...], tuple[int, ...]]):
    profile, split = entry
    complement = tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )
    reverse_profile = tuple(reversed(profile))
    return {
        entry,
        (profile, complement),
        (reverse_profile, tuple(reversed(split))),
        (reverse_profile, tuple(reversed(complement))),
    }


def find_nonzero_chain(
    plant: DirectQ4Plant, profile: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    matrices = tuple(
        plant.moments(left, right)
        for left, right in zip(profile, profile[1:])
    )
    for first in range(len(plant.supports[profile[0]])):
        seconds = np.flatnonzero(matrices[0][first])
        for second in seconds[:32]:
            thirds = np.flatnonzero(matrices[1][second])
            for third in thirds[:32]:
                fourths = np.flatnonzero(matrices[2][third])
                if fourths.size:
                    indices = first, int(second), int(third), int(fourths[0])
                    return tuple(
                        plant.supports[degree][index]
                        for degree, index in zip(profile, indices, strict=True)
                    )
    raise AssertionError(("no nonzero q4 physical chain", profile))


def q4_endpoint_and_incidence_audit(plant: DirectQ4Plant) -> tuple[int, int]:
    endpoint_parameters = set()
    for profile, split in ORBIT_REPRESENTATIVES:
        if profile == (1, 3, 5, 3):
            continue
        complement = tuple(
            degree - selected
            for degree, selected in zip(profile, split, strict=True)
        )
        endpoint_parameters.add((profile[1], split[1], split[0]))
        endpoint_parameters.add((profile[1], complement[1], complement[0]))
    endpoint_checks = 0
    for degree, selected, singleton_selected in sorted(endpoint_parameters):
        for record in (1, 3):
            direct = direct_endpoint_energy(
                plant, degree, selected, singleton_selected, record
            )
            proposed = endpoint_energy(
                degree, selected, singleton_selected, record, Q4
            )
            if direct > proposed:
                raise AssertionError(
                    (
                        "q4 endpoint energy",
                        degree,
                        selected,
                        singleton_selected,
                        record,
                        direct,
                        proposed,
                    )
                )
            endpoint_checks += 1

    signatures = tuple(
        (degree, left, right)
        for degree in (3, 5)
        for left, right in (
            (1, 1),
            (1, 3),
            (3, 1),
            (3, 3),
            (1, None),
            (3, None),
        )
    )
    incidence_checks = 0
    for degree, left, right in signatures:
        direct = direct_incidence_maxima(Q4, degree, left, right)
        proposed = block_incidence_table(degree, left, right, Q4)
        if direct != proposed:
            raise AssertionError(
                ("q4 physical incidence", degree, left, right, direct, proposed)
            )
        incidence_checks += 1
    for degree in (3, 5):
        for record in (1, 3):
            direct = direct_one_axis_incidence_maxima(Q4, degree, record)
            proposed = tuple(
                odd_record_incidence(Q4, degree, record, selected)
                for selected in range(degree + 1)
            )
            if direct != proposed:
                raise AssertionError(
                    ("q4 one-axis incidence", degree, record, direct, proposed)
                )
            incidence_checks += 1
    return endpoint_checks, incidence_checks


def q4_chain_sector_audit(plant: DirectQ4Plant) -> int:
    checks = 0
    q = Q4
    for first_record, second_record in product((1, 3), repeat=2):
        endpoint = plant.filtered_maximum(
            1, 3, (None, 1), (1, first_record)
        )
        first_link = plant.filtered_maximum(
            3,
            5,
            (1, first_record),
            (first_record, second_record),
        )
        second_link = plant.filtered_maximum(
            5,
            3,
            (first_record, second_record),
            (second_record, None),
        )
        endpoint_bound = (
            Fraction(1, q * (q - 1))
            if first_record == 1
            else Fraction(1, q)
        )
        first_bound = (
            Fraction(q + 2, q * (q - 1) * (q - 2))
            if first_record == 1
            else Fraction(3, (q - 3) * comb(q, 3))
        )
        second_bound = (
            Fraction(q + 2, q * (q - 1) * (q - 2))
            if second_record == 1
            else Fraction(3, (q - 3) * comb(q, 3))
        )
        if endpoint > endpoint_bound or first_link > first_bound or second_link > second_bound:
            raise AssertionError(
                (
                    "q4 chain link sector",
                    first_record,
                    second_record,
                    endpoint,
                    first_link,
                    second_link,
                )
            )
        for degree, left, right in (
            (3, 1, first_record),
            (5, first_record, second_record),
        ):
            physical = direct_incidence_maxima(q, degree, left, right)
            for selected, value in enumerate(physical):
                relaxed = two_axis_relaxed_incidence(
                    q, degree, left, right, selected
                )
                if value > relaxed:
                    raise AssertionError(
                        ("q4 two-axis relaxation", degree, left, right, selected)
                    )
        checks += 1
    return checks


def q4_orbit_symmetry_audit(plant: DirectQ4Plant) -> int:
    all_candidates = set(candidate_entries())
    if set().union(*(orbit(entry) for entry in ORBIT_REPRESENTATIVES)) != all_candidates:
        raise AssertionError("ten-orbit candidate reconstruction")
    checked = 0
    chains: dict[tuple[int, ...], tuple[tuple[int, ...], ...]] = {}
    for representative in ORBIT_REPRESENTATIVES:
        profile, _ = representative
        supports = chains.setdefault(profile, find_nonzero_chain(plant, profile))
        direct = Fraction(1)
        for left, right in zip(supports, supports[1:]):
            direct *= plant.moment(left, right)
        reversed_supports = tuple(
            transpose_support(support, Q4) for support in reversed(supports)
        )
        reversed_direct = Fraction(1)
        for left, right in zip(reversed_supports, reversed_supports[1:]):
            reversed_direct *= plant.moment(left, right)
        if direct != reversed_direct or not direct:
            raise AssertionError(
                ("q4 physical reversal", representative, direct, reversed_direct)
            )
        for entry in orbit(representative):
            entry_profile, split = entry
            physical_supports = (
                supports
                if entry_profile == profile
                else reversed_supports
            )
            for support, selected in zip(physical_supports, split, strict=True):
                row = set(support[:selected])
                column = set(support[selected:])
                if row.intersection(column) or row.union(column) != set(support):
                    raise AssertionError(("physical complement mask", entry))
            checked += 1
    return checked


def q8_link_and_counterexample_audit() -> tuple[int, Fraction, Fraction]:
    q = Q8
    singleton = (0,)
    cubic_l = (0, 1, q)
    cubic_row = (0, 1, 2)
    cubic_column = (0, q, 2 * q)
    quintic_zero_xor = (0, 1, 2, 3, q)

    e1 = direct_permutation_moment(q, singleton, cubic_l)
    e3 = direct_permutation_moment(q, singleton, cubic_row)
    if abs(e1) != Fraction(1, q * (q - 1)):
        raise AssertionError(("q8 endpoint record one", e1))
    if abs(e3) != Fraction(1, q):
        raise AssertionError(("q8 endpoint record three", e3))

    compatible_cubic = direct_permutation_moment(q, cubic_l, cubic_row)
    compatible_bound = Fraction(q + 2, q * (q - 1) * (q - 2))
    if abs(compatible_cubic) > compatible_bound:
        raise AssertionError(("q8 compatible cubic bound", compatible_cubic))

    generic_cubic = direct_permutation_moment(q, cubic_column, cubic_row)
    if abs(generic_cubic) != Fraction(1, q):
        raise AssertionError(("q8 generic cubic record one", generic_cubic))

    generic_cubic_quintic = direct_permutation_moment(
        q, cubic_column, quintic_zero_xor
    )
    if abs(generic_cubic_quintic) != Fraction(1, q):
        raise AssertionError(
            ("q8 generic cubic-quintic record one", generic_cubic_quintic)
        )

    # Selected representatives of every quintic endpoint row-pattern type
    # used by the exact endpoint energies: 5, 4+1 with zero/nonzero xor,
    # 3+2, and 2+2+1 with equal/distinct even-group xors.
    quintic_geometries = {
        "five": ((0, 1, 2, 3, 4), Fraction(1, q)),
        "four_one_zero": ((0, 1, 2, 3, q), Fraction(1, q)),
        "four_one_nonzero": ((0, 1, 2, 4, q), Fraction(1, q * (q - 1))),
        "three_two": ((0, 1, 2, q, q + 1), Fraction(1, q * (q - 1))),
        "two_two_one_equal": (
            (0, 1, q + 2, q + 3, 2 * q),
            Fraction(1, q * (q - 1)),
        ),
        "two_two_one_distinct": (
            (0, 1, q, q + 2, 2 * q + 3),
            Fraction(2, q * (q - 1) * (q - 2)),
        ),
    }
    for label, (support, expected_magnitude) in quintic_geometries.items():
        observed = abs(direct_permutation_moment(q, singleton, support))
        if observed != expected_magnitude:
            raise AssertionError(("q8 quintic endpoint geometry", label, observed))

    # Exact endpoint-compatible four-block counterexample in sector (1,3,1).
    cubic = cubic_row
    quintic = (0, 1, q, 2 * q, 3 * q)
    final_cubic = cubic_row
    moments = (
        direct_permutation_moment(q, singleton, cubic),
        direct_permutation_moment(q, cubic, quintic),
        direct_permutation_moment(q, quintic, final_cubic),
    )
    expected = (
        Fraction(1, q),
        -Fraction(1, (q - 3) * comb(q, 3)),
        -Fraction(1, q),
    )
    if moments != expected:
        raise AssertionError(("q8 physical counterexample moments", moments, expected))
    if tuple(support_records(support, q) for support in (cubic, quintic, final_cubic)) != (
        (1, 3),
        (3, 1),
        (1, 3),
    ):
        raise AssertionError("q8 counterexample record geometry")
    actual = abs(moments[0] * moments[1] * moments[2])
    m1 = Fraction(q + 2, q * (q - 1) * (q - 2))
    m3 = Fraction(3, (q - 3) * comb(q, 3))
    claimed = Fraction(1, q) * m3 * m1
    if not actual > claimed or (actual, claimed) != (
        Fraction(1, 17920),
        Fraction(1, 25088),
    ):
        raise AssertionError(("q8 claimed maximum survives", actual, claimed))

    # The same explicitly summed permutation family invalidates the q=64
    # maximum by a factor 217/11.  All quantities remain exact Fractions.
    q64 = 64
    actual_q64 = Fraction(1, q64**2 * (q64 - 3) * comb(q64, 3))
    claimed_q64 = (
        Fraction(1, q64)
        * Fraction(3, (q64 - 3) * comb(q64, 3))
        * Fraction(q64 + 2, q64 * (q64 - 1) * (q64 - 2))
    )
    if actual_q64 / claimed_q64 != Fraction(217, 11):
        raise AssertionError(("q64 counterexample ratio", actual_q64, claimed_q64))
    return 13, actual, claimed


def endpoint_link_scope_audit() -> int:
    q = 64
    # Profile (1,3,3,5): the first cubic has incoming record one, so the
    # refined cubic--cubic record-one bound is in scope.
    if endpoint_remaining_link_bound((1, 3, 3, 5), (1, 1, 1), 1, q) != Fraction(
        q + 2, q * (q - 1) * (q - 2)
    ):
        raise AssertionError("lost compatible cubic refinement")
    # Profile (1,5,3,3), sector (1,3,1): the left cubic's incoming record is
    # three, so the exact vertical/horizontal example forces the generic 1/q.
    if endpoint_remaining_link_bound((1, 5, 3, 3), (1, 3, 1), 2, q) != Fraction(
        1, q
    ):
        raise AssertionError("unsafe cubic refinement retained")
    return 2


def numerical_artifact_audit() -> tuple[int, float]:
    artifact = loads(
        (
            ROOT
            / "artifacts"
            / "q64_masked_recovered_cubic_quintic_incidence_repair.json"
        ).read_text(encoding="utf-8")
    )
    if artifact["schema"] != "round4_q64_masked_recovered_cubic_quintic_incidence_repair_v3":
        raise AssertionError("stale recovered audit artifact schema")
    rows = artifact["repaired_registry_entries"]
    if len(rows) != 28 or len(artifact["rejected_chain_entries"]) != 12:
        raise AssertionError("recovered artifact inventory")
    getcontext().prec = 100
    maximum = 0.0
    for row in rows:
        exact_sum = Decimal(0)
        for sector in row["sector_squared_coefficients"]:
            exact = Decimal(sector["numerator"]) / Decimal(sector["denominator"])
            exact_sum += exact.sqrt()
        displayed = Decimal.from_float(row["outward_coefficient"])
        if displayed < exact_sum:
            raise AssertionError(
                ("artifact coefficient rounded inward", row["profile"], row["split"])
            )
        maximum = max(maximum, row["outward_coefficient"])
    if maximum >= 1 or 1 - maximum < 1e-3:
        raise AssertionError(("insufficient coefficient margin", maximum))
    counterexample = artifact["counterexample"]
    if Fraction(
        counterexample["exact_entry_numerator"],
        counterexample["exact_entry_denominator"],
    ) <= Fraction(
        counterexample["claimed_maximum_entry_numerator"],
        counterexample["claimed_maximum_entry_denominator"],
    ):
        raise AssertionError("artifact lost rejecting counterexample")
    return len(rows), maximum


def registry_audit() -> tuple[int, int]:
    candidates = set(candidate_entries())
    repaired = set(repaired_entries())
    rejected = set(rejected_chain_entries())
    if len(candidates) != 40 or len(repaired) != 28 or len(rejected) != 12:
        raise AssertionError("recovered independent registry counts")
    if repaired.intersection(rejected) or repaired.union(rejected) != candidates:
        raise AssertionError("recovered independent registry partition")
    if any(canonical_entry(entry)[0] != (1, 3, 5, 3) for entry in rejected):
        raise AssertionError("wrong rejected recovered family")
    return len(repaired), len(rejected)


def main() -> None:
    plant = DirectQ4Plant()
    endpoint_checks, incidence_checks = q4_endpoint_and_incidence_audit(plant)
    sector_checks = q4_chain_sector_audit(plant)
    orbit_checks = q4_orbit_symmetry_audit(plant)
    q8_checks, actual, claimed = q8_link_and_counterexample_audit()
    scope_checks = endpoint_link_scope_audit()
    artifact_rows, maximum = numerical_artifact_audit()
    repaired, rejected = registry_audit()
    print(
        "q64 recovered cubic-quintic independent audit passed: "
        "verdict=rejected_complete_40_claim,"
        f"retained={repaired},quarantined={rejected},"
        f"q4_endpoint={endpoint_checks},q4_incidence={incidence_checks},"
        f"q4_sectors={sector_checks},q4_orientations={orbit_checks},"
        f"q8_links={q8_checks},scope={scope_checks},artifact_rows={artifact_rows},"
        f"counterexample={actual}>{claimed},maximum_retained={maximum:.12g}"
    )


if __name__ == "__main__":
    main()
