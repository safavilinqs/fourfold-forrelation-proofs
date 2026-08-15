#!/usr/bin/env python3
"""Independent physical audit of the final eighty q64 residual entries."""

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

import q64_final_residual_chain_contraction as theorem  # noqa: E402
from q64_degree_ten_completion_row_insertion import orbit  # noqa: E402
from q64_masked_local_walsh_repair import mechanism as local_mechanism  # noqa: E402
from q64_same_side_whole_link_insertion import remaining_entries  # noqa: E402


Q4 = 4
Q8 = 8
CHAIN_GENERATORS = (
    ((1, 3, 3, 5), (0, 3, 0, 3)),
    ((1, 3, 5, 1), (0, 3, 1, 1)),
    ((1, 3, 5, 1), (0, 3, 2, 0)),
    ((1, 3, 7, 1), (0, 3, 2, 1)),
    ((1, 3, 7, 1), (0, 3, 3, 0)),
    ((1, 3, 5, 3), (0, 3, 3, 0)),
    ((1, 5, 3, 3), (0, 3, 0, 3)),
    ((1, 5, 3, 3), (0, 3, 3, 0)),
)


def walsh(row: int, column: int) -> int:
    return -1 if (row & column).bit_count() % 2 else 1


def support_records(support: tuple[int, ...], order: int) -> tuple[int, int]:
    records = []
    for axis in (0, 1):
        counts = Counter(divmod(coordinate, order)[axis] for coordinate in support)
        records.append(sum(value % 2 for value in counts.values()))
    return tuple(records)  # type: ignore[return-value]


def direct_permutation_moment(
    order: int,
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> Fraction:
    """Average one link directly over permutations after the exact sign sum."""

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
    """Exact q=4 signed-permutation plant, independent of moment formulas.

    The largest stored feature array is 11,440 by 384 int16 entries, below
    9 MB.  Dense matrices are limited to at most three million entries; the
    degree-seven links are evaluated one physical pair at a time.
    """

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
            for degree in (1, 3, 5, 7)
        }
        self.support_index = {
            degree: {support: index for index, support in enumerate(supports)}
            for degree, supports in self.supports.items()
        }
        self.indices_by_record = {
            (degree, records): tuple(
                index
                for index, support in enumerate(self.supports[degree])
                if support_records(support, q) == records
            )
            for degree in self.supports
            for records in product(range(1, degree + 1, 2), repeat=2)
        }
        self._features: dict[tuple[str, int], np.ndarray] = {}
        self._matrices: dict[tuple[int, int], np.ndarray] = {}

    def features(self, side: str, degree: int) -> np.ndarray:
        key = side, degree
        if key not in self._features:
            values = self.left_values if side == "left" else self.right_values
            indices = np.asarray(self.supports[degree], dtype=np.int16)
            self._features[key] = np.prod(
                values[:, indices], axis=2, dtype=np.int16
            ).T
        return self._features[key]

    def matrix(self, left_degree: int, right_degree: int) -> np.ndarray:
        key = left_degree, right_degree
        entries = len(self.supports[left_degree]) * len(self.supports[right_degree])
        if entries > 3_000_000:
            raise ValueError(("dense q4 link deliberately disabled", key, entries))
        if key not in self._matrices:
            self._matrices[key] = (
                self.features("left", left_degree)
                @ self.features("right", right_degree).T
            )
        return self._matrices[key]

    def numerator(self, left: tuple[int, ...], right: tuple[int, ...]) -> int:
        left_index = self.support_index[len(left)][left]
        right_index = self.support_index[len(right)][right]
        if len(self.supports[len(left)]) * len(self.supports[len(right)]) <= 3_000_000:
            return int(self.matrix(len(left), len(right))[left_index, right_index])
        return int(
            self.features("left", len(left))[left_index]
            @ self.features("right", len(right))[right_index]
        )

    def endpoint_energy(self, support: tuple[int, ...], endpoint_first: bool) -> Fraction:
        total = 0
        for singleton in self.supports[1]:
            numerator = (
                self.numerator(singleton, support)
                if endpoint_first
                else self.numerator(support, singleton)
            )
            total += numerator * numerator
        return Fraction(total, self.group_size**2)

    def completions(
        self,
        degree: int,
        partial: tuple[int, ...],
        row_record: int | None = None,
        column_record: int | None = None,
    ) -> tuple[tuple[int, ...], ...]:
        fixed = set(partial)
        return tuple(
            support
            for support in self.supports[degree]
            if fixed.issubset(support)
            and (
                row_record is None
                or support_records(support, Q4)[0] == row_record
            )
            and (
                column_record is None
                or support_records(support, Q4)[1] == column_record
            )
        )


def direct_incidence_maximum(
    plant: DirectQ4Plant,
    degree: int,
    selected: int,
    row_record: int | None,
    column_record: int | None,
) -> int:
    counts: Counter[tuple[int, ...]] = Counter()
    for support in plant.supports[degree]:
        observed = support_records(support, Q4)
        if row_record is not None and observed[0] != row_record:
            continue
        if column_record is not None and observed[1] != column_record:
            continue
        counts.update(combinations(support, selected))
    return max(counts.values(), default=0)


def direct_fixed_endpoint_partial_energy(
    plant: DirectQ4Plant,
    degree: int,
    selected: int,
    row_record: int,
) -> Fraction:
    """Maximum fixed-endpoint energy through a fixed physical partial."""

    accumulators: defaultdict[tuple[tuple[int, ...], tuple[int, ...]], int]
    accumulators = defaultdict(int)
    for support in plant.supports[degree]:
        if support_records(support, Q4)[0] != row_record:
            continue
        for endpoint in plant.supports[1]:
            numerator = plant.numerator(endpoint, support)
            if not numerator:
                continue
            for partial in combinations(support, selected):
                accumulators[endpoint, partial] += numerator * numerator
    return Fraction(max(accumulators.values(), default=0), plant.group_size**2)


def direct_variable_endpoint_partial_energy(
    plant: DirectQ4Plant,
    degree: int,
    selected: int,
    row_record: int,
    column_record: int,
) -> Fraction:
    """Maximum variable terminal-endpoint energy through a fixed partial."""

    accumulators: defaultdict[tuple[int, ...], int] = defaultdict(int)
    for support in plant.supports[degree]:
        if support_records(support, Q4) != (row_record, column_record):
            continue
        energy_numerator = sum(
            plant.numerator(support, endpoint) ** 2
            for endpoint in plant.supports[1]
        )
        for partial in combinations(support, selected):
            accumulators[partial] += energy_numerator
    return Fraction(max(accumulators.values(), default=0), plant.group_size**2)


def representative_support(
    plant: DirectQ4Plant,
    degree: int,
    records: tuple[int, int],
) -> tuple[int, ...]:
    indices = plant.indices_by_record[degree, records]
    if not indices:
        raise AssertionError(("empty q4 physical record family", degree, records))
    return plant.supports[degree][indices[0]]


def q4_terminal_high_audit(plant: DirectQ4Plant) -> int:
    """Direct complete rows for singleton--cubic--cubic--quintic."""

    q = Q4
    group = plant.group_size
    b = Fraction(q + 2, q * (q - 1) * (q - 2))
    g = Fraction(1, comb(q, 3))
    checks = 0
    for first_record, second_record in product((1, 3), repeat=2):
        cubic = representative_support(plant, 3, (1, first_record))
        quintic = representative_support(plant, 5, (second_record, 1))
        partial = quintic[:3]
        endpoint_numerator = sum(
            plant.numerator(singleton, cubic) ** 2
            for singleton in plant.supports[1]
        )
        chain_numerator = 0
        for middle_index in plant.indices_by_record[3, (first_record, second_record)]:
            middle = plant.supports[3][middle_index]
            left_numerator = plant.numerator(cubic, middle)
            if not left_numerator:
                continue
            for high in plant.completions(5, partial, second_record, None):
                right_numerator = plant.numerator(middle, high)
                chain_numerator += left_numerator**2 * right_numerator**2
        direct = Fraction(endpoint_numerator * chain_numerator, group**6)

        endpoint_max = max(
            plant.endpoint_energy(
                plant.supports[3][index], endpoint_first=True
            )
            for index in plant.indices_by_record[3, (1, first_record)]
        )
        cubic_count = len(plant.indices_by_record[3, (first_record, second_record)])
        quintic_count = direct_incidence_maximum(
            plant, 5, 3, second_record, None
        )
        if (first_record, second_record) != (3, 1):
            first_bound = b if first_record == 1 else g
            second_bound = b if second_record == 1 else g
            proposed = endpoint_max * cubic_count * quintic_count * (
                first_bound * second_bound
            ) ** 2
        else:
            no_even = sum(
                len({divmod(cell, q)[1] for cell in plant.supports[3][index]}) == 1
                for index in plant.indices_by_record[3, (3, 1)]
            )
            active = cubic_count - no_even
            fixed_endpoint = direct_fixed_endpoint_partial_energy(
                plant, 5, 3, 1
            )
            proposed = endpoint_max * g**2 * (
                no_even * fixed_endpoint
                + active * quintic_count * b**2
            )
        if direct > proposed:
            raise AssertionError(
                (
                    "q4 terminal-high complete row",
                    first_record,
                    second_record,
                    direct,
                    proposed,
                )
            )
        checks += 1
    return checks


def q4_double_endpoint_audit(plant: DirectQ4Plant) -> int:
    """Direct rows for all four cubic--higher--singleton orbit types."""

    q = Q4
    group = plant.group_size
    b = Fraction(q + 2, q * (q - 1) * (q - 2))
    g = Fraction(1, comb(q, 3))
    parameters = ((5, 1, 1), (5, 2, 0), (7, 2, 1), (7, 3, 0))
    checks = 0
    for degree, selected, endpoint_selected in parameters:
        for record in (1, 3):
            cubic = representative_support(plant, 3, (1, record))
            high = representative_support(plant, degree, (record, 1))
            partial = high[:selected]
            terminal = plant.supports[1][0]
            endpoint_numerator = sum(
                plant.numerator(singleton, cubic) ** 2
                for singleton in plant.supports[1]
            )
            remainder_numerator = 0
            for completion in plant.completions(
                degree, partial, record, 1
            ):
                middle = plant.numerator(cubic, completion)
                if endpoint_selected:
                    terminal_energy = plant.numerator(completion, terminal) ** 2
                else:
                    terminal_energy = sum(
                        plant.numerator(completion, singleton) ** 2
                        for singleton in plant.supports[1]
                    )
                remainder_numerator += middle**2 * terminal_energy
            direct = Fraction(endpoint_numerator * remainder_numerator, group**6)

            endpoint_max = max(
                plant.endpoint_energy(
                    plant.supports[3][index], endpoint_first=True
                )
                for index in plant.indices_by_record[3, (1, record)]
            )
            middle_bound = (
                (b if record == 1 else g)
                if degree == 5
                else (Fraction(1, q) if record == 1 else g)
            )
            if endpoint_selected:
                terminal_bound = Fraction(
                    direct_incidence_maximum(
                        plant, degree, selected, record, 1
                    ),
                    q * q,
                )
            elif degree == 5:
                terminal_bound = Fraction(
                    direct_incidence_maximum(
                        plant, degree, selected, record, 1
                    )
                )
            elif record == 1:
                terminal_bound = direct_variable_endpoint_partial_energy(
                    plant, degree, selected, record, 1
                )
            else:
                direct_energy = direct_variable_endpoint_partial_energy(
                    plant, degree, selected, record, 1
                )
                # The q>=8 active-character inequality used by the theorem
                # deliberately excludes q=4: its three-active-group factor
                # is 3/((q-1)(q-3)), which is not at most 1/(q-1) at q=4.
                # Use the directly enumerated energy here to isolate and
                # validate the complete-row norm implication.  The q=8
                # exceptional-family test below validates the theorem's
                # actual endpoint threshold.
                terminal_bound = direct_energy
            proposed = endpoint_max * middle_bound**2 * terminal_bound
            if direct > proposed:
                raise AssertionError(
                    (
                        "q4 double-endpoint complete row",
                        degree,
                        selected,
                        endpoint_selected,
                        record,
                        direct,
                        proposed,
                    )
                )
            checks += 1
    return checks


def q4_endpoint_cubic_high_audit(plant: DirectQ4Plant) -> int:
    q = Q4
    group = plant.group_size
    b = Fraction(q + 2, q * (q - 1) * (q - 2))
    g = Fraction(1, comb(q, 3))
    checks = 0
    for record in (1, 3):
        cubic = representative_support(plant, 3, (1, record))
        high = representative_support(plant, 5, (record, 1))
        partial = high[:3]
        endpoint_numerator = sum(
            plant.numerator(singleton, cubic) ** 2
            for singleton in plant.supports[1]
        )
        middle_numerator = sum(
            plant.numerator(cubic, completion) ** 2
            for completion in plant.completions(5, partial, record, None)
        )
        direct = Fraction(endpoint_numerator * middle_numerator, group**4)
        endpoint_max = max(
            plant.endpoint_energy(
                plant.supports[3][index], endpoint_first=True
            )
            for index in plant.indices_by_record[3, (1, record)]
        )
        incidence = direct_incidence_maximum(plant, 5, 3, record, None)
        proposed = endpoint_max * incidence * (b if record == 1 else g) ** 2
        if direct > proposed:
            raise AssertionError(
                ("q4 endpoint-cubic-high row", record, direct, proposed)
            )
        checks += 1
    return checks


def q4_endpoint_high_cubic_audit(plant: DirectQ4Plant) -> int:
    """Direct checks of both complementary endpoint--quintic orientations."""

    q = Q4
    group = plant.group_size
    checks = 0
    endpoint = plant.supports[1][0]
    for record in (1, 3):
        cubic = representative_support(plant, 3, (record, 1))
        high = representative_support(plant, 5, (1, record))
        partial = high[:2]
        numerator = sum(
            plant.numerator(endpoint, completion) ** 2
            * plant.numerator(completion, cubic) ** 2
            for completion in plant.completions(5, partial, 1, record)
        )
        direct = Fraction(numerator, group**4)
        incidence = direct_incidence_maximum(plant, 5, 2, 1, record)
        proposed = Fraction(incidence, q * q * comb(q, record) ** 2)
        if direct > proposed:
            raise AssertionError(
                ("q4 endpoint-high-cubic completion", record, direct, proposed)
            )
        checks += 1

    endpoint_slice = direct_fixed_endpoint_partial_energy(plant, 5, 3, 1)
    if endpoint_slice != Fraction(7, 8):
        raise AssertionError(("q4 fixed-endpoint quintic triple slice", endpoint_slice))
    for cubic_record in (1, 3):
        # The endpoint label is a normalized Walsh Schur dressing.  It must
        # not be summed as q^2 unrelated completion columns; the scalar
        # feature is the fixed-endpoint slice checked above.  Enumerate the
        # remaining physical quintic--cubic factor separately.  At q=4 the
        # record-three correction is larger than its q>=8 theorem value, so
        # the test records the exact small-order maximum rather than applying
        # an out-of-range formula.
        high_indices = [
            index
            for index, support in enumerate(plant.supports[5])
            if support_records(support, q)[0] == 1
            and support_records(support, q)[1] == cubic_record
        ]
        cubic_indices = plant.indices_by_record[3, (cubic_record, 1)]
        link_numerator = int(
            np.abs(plant.matrix(5, 3)[np.ix_(high_indices, cubic_indices)]).max(
                initial=0
            )
        )
        link_bound = Fraction(link_numerator, group)
        if not 0 < link_bound <= 1:
            raise AssertionError(("q4 endpoint-high-cubic dressing", cubic_record, link_bound))
        checks += 1
    return checks


def q8_link_geometry_audit() -> int:
    """Direct permutation sums for every exceptional link geometry."""

    q = Q8
    b = Fraction(q + 2, q * (q - 1) * (q - 2))
    g = Fraction(1, comb(q, 3))
    endpoint_cubic = tuple(sorted((0 * q + 0, 0 * q + 1, 1 * q + 0)))
    active_cubic = tuple(sorted((0 * q + 0, 1 * q + 0, 2 * q + 1)))
    no_even_cubic = tuple(sorted((0 * q + 0, 1 * q + 0, 2 * q + 0)))
    quintics = {
        "221": tuple(sorted((0 * q + 0, 0 * q + 1, 1 * q + 0, 1 * q + 2, 2 * q + 3))),
        "32": tuple(sorted((0 * q + 0, 0 * q + 1, 0 * q + 2, 1 * q + 0, 1 * q + 1))),
        "41_zero": tuple(sorted((0 * q + 0, 0 * q + 1, 0 * q + 2, 0 * q + 3, 1 * q + 0))),
        "41_active": tuple(sorted((0 * q + 0, 0 * q + 1, 0 * q + 2, 0 * q + 4, 1 * q + 0))),
    }
    endpoint = direct_permutation_moment(q, (0,), endpoint_cubic)
    if abs(endpoint) != Fraction(1, q * (q - 1)):
        raise AssertionError(("q8 endpoint-compatible cubic", endpoint))
    checks = 1
    for label, quintic in quintics.items():
        active = abs(direct_permutation_moment(q, active_cubic, quintic))
        if active > b:
            raise AssertionError(("q8 active record-one link", label, active, b))
        no_even = abs(direct_permutation_moment(q, no_even_cubic, quintic))
        if no_even > Fraction(1, q):
            raise AssertionError(("q8 no-even endpoint reduction", label, no_even))
        checks += 2

    record_three_cubic = tuple(sorted((0, 1, 2)))
    record_three_quintic = tuple(sorted((0, 1, q, 2 * q, 3 * q)))
    observed = abs(
        direct_permutation_moment(q, record_three_cubic, record_three_quintic)
    )
    if observed > g:
        raise AssertionError(("q8 record-three universal link", observed, g))
    return checks + 1


def q8_zero_active_septimic_audit() -> tuple[int, int]:
    """Enumerate the exceptional 4+3 record-(3,1) septimic family."""

    q = Q8
    zero_four_sets = tuple(
        support
        for support in combinations(range(q), 4)
        if support[0] ^ support[1] ^ support[2] ^ support[3] == 0
    )
    supports = set()
    for degree_four_column in range(q):
        for degree_three_column in range(q):
            if degree_four_column == degree_three_column:
                continue
            for four_rows in zero_four_sets:
                four_set = set(four_rows)
                for shared_rows in combinations(four_rows, 2):
                    for outside in set(range(q)) - four_set:
                        three_rows = (*shared_rows, outside)
                        support = tuple(
                            sorted(
                                tuple(row * q + degree_four_column for row in four_rows)
                                + tuple(row * q + degree_three_column for row in three_rows)
                            )
                        )
                        supports.add(support)
    expected_total = q * (q - 1) * (q * (q - 1) * (q - 2) // 24) * 6 * (q - 4)
    if len(supports) != expected_total:
        raise AssertionError(("q8 zero-active total", len(supports), expected_total))
    if any(support_records(support, q) != (3, 1) for support in supports):
        raise AssertionError("q8 zero-active record classification")
    cell_counts: Counter[int] = Counter()
    triple_counts: Counter[tuple[int, ...]] = Counter()
    for support in supports:
        cell_counts.update(support)
        triple_counts.update(combinations(support, 3))
    expected_cell = 7 * expected_total // (q * q)
    if set(cell_counts.values()) != {expected_cell}:
        raise AssertionError(("q8 zero-active transitivity", set(cell_counts.values())))
    if max(triple_counts.values()) > expected_cell:
        raise AssertionError("q8 zero-active fixed-triple domination")

    zero_support = tuple(
        sorted(
            tuple(row * q for row in (0, 1, 2, 3))
            + tuple(row * q + 1 for row in (0, 1, 4))
        )
    )
    zero_energy = sum(
        direct_permutation_moment(q, zero_support, (endpoint,)) ** 2
        for endpoint in range(q * q)
    )
    if zero_energy != 1:
        raise AssertionError(("q8 zero-active endpoint energy", zero_energy))

    active_support = tuple(sorted((0, 8, 16, 1, 25, 10, 34)))
    if support_records(active_support, q) != (3, 1):
        raise AssertionError("q8 active septimic records")
    active_energy = sum(
        direct_permutation_moment(q, active_support, (endpoint,)) ** 2
        for endpoint in range(q * q)
    )
    if active_energy > Fraction(1, (q - 1) ** 2):
        raise AssertionError(("q8 active septimic endpoint energy", active_energy))
    return expected_total, expected_cell


def inventory_and_numerical_audit() -> tuple[int, int, float]:
    entries = set(remaining_entries())
    local = {entry for entry in entries if local_mechanism(entry) is not None}
    chain = entries - local
    expected_chain = set().union(*(orbit(entry) for entry in CHAIN_GENERATORS))
    if len(entries) != 80 or len(local) != 48 or len(chain) != 32:
        raise AssertionError(("final residual inventory", len(entries), len(local), len(chain)))
    if chain != expected_chain or len({min(orbit(entry)) for entry in chain}) != 8:
        raise AssertionError("eight-orbit chain reconstruction")

    artifact_path = ROOT / "artifacts" / "q64_final_residual_chain_contraction.json"
    committed = artifact_path.read_text(encoding="utf-8")
    artifact = loads(committed)
    if artifact["schema"] != "round4_q64_final_residual_chain_contraction_v1":
        raise AssertionError("stale final-residual artifact schema")
    rows = artifact["registry_entries"]
    if len(rows) != 80:
        raise AssertionError(("final-residual artifact inventory", len(rows)))
    getcontext().prec = 100
    maximum = 0.0
    for row in rows:
        exact = Decimal(row["squared_coefficient_numerator"]) / Decimal(
            row["squared_coefficient_denominator"]
        )
        displayed = Decimal.from_float(row["outward_coefficient"])
        if displayed * displayed < exact:
            raise AssertionError(
                ("final-residual coefficient rounded inward", row["profile"], row["split"])
            )
        maximum = max(maximum, row["outward_coefficient"])
    if maximum >= 1 or 1 - maximum < 0.5:
        raise AssertionError(("final-residual coefficient margin", maximum))
    if committed != theorem.artifact_text(theorem.diagnostic()):
        raise AssertionError("stale final-residual theorem artifact")
    return len(local), len(chain), maximum


def main() -> None:
    plant = DirectQ4Plant()
    terminal_checks = q4_terminal_high_audit(plant)
    double_endpoint_checks = q4_double_endpoint_audit(plant)
    endpoint_cubic_checks = q4_endpoint_cubic_high_audit(plant)
    endpoint_high_checks = q4_endpoint_high_cubic_audit(plant)
    q8_link_checks = q8_link_geometry_audit()
    zero_total, zero_cell = q8_zero_active_septimic_audit()
    local, chain, maximum = inventory_and_numerical_audit()
    print(
        "q64 final residual independent audit passed: "
        f"inventory={local}+{chain},q4_terminal={terminal_checks},"
        f"q4_double_endpoint={double_endpoint_checks},"
        f"q4_endpoint_cubic={endpoint_cubic_checks},"
        f"q4_endpoint_high={endpoint_high_checks},q8_links={q8_link_checks},"
        f"q8_zero_family={zero_total}/{zero_cell},maximum={maximum:.12g}"
    )


if __name__ == "__main__":
    main()
