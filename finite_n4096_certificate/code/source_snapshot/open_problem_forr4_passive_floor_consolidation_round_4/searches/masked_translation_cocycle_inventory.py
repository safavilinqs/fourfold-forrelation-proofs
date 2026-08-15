#!/usr/bin/env python3
"""Classify the projective translation cocycles of the 97 templates."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from json import dumps
from math import log2
from pathlib import Path

from masked_translation_covariance import translation_cocycle
from q4_masked_template_screen import template_inventory


ROOT = Path(__file__).resolve().parents[1]
Shift = tuple[int, int]
FourShifts = tuple[Shift, Shift, Shift, Shift]


@dataclass(frozen=True)
class CocycleTemplateRow:
    class_name: str
    profile: tuple[int, ...]
    split: tuple[int, ...]
    split_parities: tuple[int, ...]
    row_commutator_rank: int
    column_commutator_rank: int
    row_column_cocycles_equal: bool
    normalized_rank: int


def predicted_normalized_rank(parities: tuple[int, ...]) -> int:
    """Closed form for the four-block projective commutator rank."""

    if parities in ((0, 0, 0, 0), (1, 1, 1, 1)):
        return 0
    if parities in ((0, 0, 1, 1), (1, 1, 0, 0)):
        return 4
    if sum(parities) == 2:
        return 8
    raise ValueError(("unexpected split parity pattern", parities))


@dataclass(frozen=True)
class CocycleInventory:
    order: int
    axis_bits: int
    translation_bits: int
    templates: int
    rank_counts: dict[int, int]
    normalized_rank_counts: dict[int, int]
    parity_rank_counts: dict[str, int]
    rows: tuple[CocycleTemplateRow, ...]


def basis_shift(order: int, bit: int) -> FourShifts:
    axis_bits = int(log2(order))
    coordinate_bits = 2 * axis_bits
    value = 1 << bit
    result = []
    for block in range(4):
        block_value = (value >> (block * coordinate_bits)) & (
            (1 << coordinate_bits) - 1
        )
        result.append(
            (block_value >> axis_bits, block_value & (order - 1))
        )
    return tuple(result)


def gf2_rank(rows: list[int], dimension: int) -> int:
    rank = 0
    for column in range(dimension):
        pivot = next(
            (
                candidate
                for candidate in range(rank, dimension)
                if (rows[candidate] >> column) & 1
            ),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for candidate in range(dimension):
            if candidate != rank and (rows[candidate] >> column) & 1:
                rows[candidate] ^= rows[rank]
        rank += 1
    return rank


def commutator_rank(order: int, block_sizes: tuple[int, ...]) -> int:
    axis_bits = int(log2(order))
    dimension = 8 * axis_bits
    basis = tuple(basis_shift(order, bit) for bit in range(dimension))
    rows = []
    for left in basis:
        row = 0
        for column, right in enumerate(basis):
            commutator = translation_cocycle(
                order, block_sizes, left, right
            ) * translation_cocycle(order, block_sizes, right, left)
            if commutator < 0:
                row |= 1 << column
        rows.append(row)
    return gf2_rank(rows, dimension)


def cocycle_basis_table(
    order: int, block_sizes: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    axis_bits = int(log2(order))
    dimension = 8 * axis_bits
    basis = tuple(basis_shift(order, bit) for bit in range(dimension))
    return tuple(
        tuple(
            translation_cocycle(order, block_sizes, left, right)
            for right in basis
        )
        for left in basis
    )


def diagnostic(order: int) -> CocycleInventory:
    axis_bits = int(log2(order))
    rows = []
    for class_name, (profile, split) in template_inventory():
        column = tuple(
            degree - selected
            for degree, selected in zip(profile, split, strict=True)
        )
        row_rank = commutator_rank(order, split)
        column_rank = commutator_rank(order, column)
        if row_rank != column_rank:
            raise AssertionError(
                ("row/column projective rank mismatch", profile, split)
            )
        cocycles_equal = cocycle_basis_table(
            order, split
        ) == cocycle_basis_table(order, column)
        if not cocycles_equal:
            raise AssertionError(
                ("row/column projective cocycle mismatch", profile, split)
            )
        if row_rank % axis_bits:
            raise AssertionError(("non-scaling cocycle rank", order, split))
        parities = tuple(value % 2 for value in split)
        normalized_rank = row_rank // axis_bits
        predicted_rank = predicted_normalized_rank(parities)
        if normalized_rank != predicted_rank:
            raise AssertionError(
                (
                    "cocycle rank disagrees with parity formula",
                    order,
                    split,
                    normalized_rank,
                    predicted_rank,
                )
            )
        rows.append(
            CocycleTemplateRow(
                class_name=class_name,
                profile=profile,
                split=split,
                split_parities=parities,
                row_commutator_rank=row_rank,
                column_commutator_rank=column_rank,
                row_column_cocycles_equal=cocycles_equal,
                normalized_rank=normalized_rank,
            )
        )
    rank_counts = Counter(row.row_commutator_rank for row in rows)
    normalized_counts = Counter(row.normalized_rank for row in rows)
    parity_counts = Counter(
        f"{''.join(map(str, row.split_parities))}:"
        f"{row.normalized_rank}"
        for row in rows
    )
    return CocycleInventory(
        order=order,
        axis_bits=axis_bits,
        translation_bits=8 * axis_bits,
        templates=len(rows),
        rank_counts=dict(sorted(rank_counts.items())),
        normalized_rank_counts=dict(sorted(normalized_counts.items())),
        parity_rank_counts=dict(sorted(parity_counts.items())),
        rows=tuple(rows),
    )


def artifact_text(result: CocycleInventory) -> str:
    payload = {
        "schema": "round4_masked_translation_cocycle_inventory_v1",
        "result": asdict(result),
        "evidence_label": (
            "exact GF(2) commutator ranks of the projective masked-translation "
            "action; row and complementary-column ranks agree for every one "
            "of the 97 structural templates"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, choices=(4, 8, 64), default=4)
    parser.add_argument("--write-artifact", action="store_true")
    arguments = parser.parse_args()
    result = diagnostic(arguments.order)
    if arguments.write_artifact:
        path = (
            ROOT
            / "artifacts"
            / f"q{arguments.order}_masked_translation_cocycle_inventory.json"
        )
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "masked translation cocycle inventory: "
        f"q={result.order},"
        f"templates={result.templates},"
        f"ranks={result.rank_counts},"
        f"normalized_ranks={result.normalized_rank_counts},"
        "status=three_projective_clifford_types"
    )


if __name__ == "__main__":
    main()
