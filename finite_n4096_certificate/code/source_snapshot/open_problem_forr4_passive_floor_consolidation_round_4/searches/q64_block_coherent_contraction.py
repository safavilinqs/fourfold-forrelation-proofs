#!/usr/bin/env python3
"""Close every q64 block-coherent balanced high-sector entry.

For a whole-block occurrence cut, the accepted weighted three-link path
contraction applies without any within-block disjointness masks.  Decompose
the signed-permutation moment into compatible odd-record triples.  A link
with record size r has maximum entry 1/C(q,r), while its pure-record operator
norm is at most q^r/(q)_r.  The complete cut table then gives one exact
rational coefficient per record triple.  Summing the sectors closes all 70
block-coherent entries in the open degree-ten/twelve profiles at q=64.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from json import dumps
from math import comb, prod
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
ROUND3_SEARCHES = (
    ROOT.parent / "open_problem_forr4_passive_floor_consolidation_round_3" / "searches"
)
sys.path.insert(0, str(ROUND3_SEARCHES))

from high_degree_record_incidence_frontier import (  # noqa: E402
    compatible_record_triples,
)

from q64_paper_target_gate import (  # noqa: E402
    DIMENSION,
    MODES,
    ORDER,
    OptimizedLedger,
    balanced_open_entries,
    optimize,
    two_tier_coefficients,
)
import occupation_compatible_sector_optimization as occupation  # noqa: E402


Profile = tuple[int, ...]
Split = tuple[int, ...]
Records = tuple[int, int, int]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class BlockCoherentResult:
    order: int
    dimension: int
    sign_modes: int
    closed_entries: int
    record_sector_bounds: int
    minimum_coefficient: float
    minimum_coefficient_exact: str
    maximum_coefficient: float
    maximum_coefficient_exact: str
    maximum_entry_profile: Profile
    maximum_entry_split: Split
    previous_two_tier: OptimizedLedger
    block_coherent_inserted: OptimizedLedger
    margin_improvement: float


def falling_factorial(order: int, size: int) -> int:
    return prod(range(order - size + 1, order + 1))


def canonical_cut(profile: Profile, split: Split) -> frozenset[int]:
    if not all(
        selected in (0, degree) for degree, selected in zip(profile, split, strict=True)
    ):
        raise ValueError(("not block coherent", profile, split))
    selected = frozenset(
        block
        for block, (degree, count) in enumerate(zip(profile, split, strict=True))
        if count == degree
    )
    complement = frozenset(range(4)) - selected
    return min(selected, complement, key=lambda value: (len(value), sorted(value)))


def sector_bound(
    profile: Profile, split: Split, records: Records, order: int = ORDER
) -> Fraction:
    """Profile-aware wrapper for the record-sector cut table."""

    kappa = tuple(Fraction(1, comb(order, record)) for record in records)
    operator = tuple(
        Fraction(order**record, falling_factorial(order, record)) for record in records
    )
    mask = canonical_cut(profile, split)
    if not mask:
        return prod(kappa, start=Fraction(1))
    if mask == frozenset({0}):
        return kappa[1] * kappa[2]
    if mask == frozenset({1}):
        return kappa[2]
    if mask == frozenset({2}):
        return kappa[0]
    if mask == frozenset({3}):
        return kappa[0] * kappa[1]
    if mask in (frozenset({0, 1}), frozenset({2, 3})):
        return kappa[0] * kappa[2]
    if mask in (frozenset({1, 2}), frozenset({0, 3})):
        return kappa[1]
    if mask in (frozenset({0, 2}), frozenset({1, 3})):
        return operator[0] * kappa[1] * operator[2]
    raise AssertionError(("unclassified block cut", mask))


def block_coherent_entries() -> tuple[ProfileSplit, ...]:
    return tuple(
        (profile, split)
        for profile, split in balanced_open_entries()
        if all(
            selected in (0, degree)
            for degree, selected in zip(profile, split, strict=True)
        )
    )


def block_coherent_coefficient(
    profile: Profile, split: Split, order: int = ORDER
) -> Fraction:
    return sum(
        (
            sector_bound(profile, split, records, order)
            for records in compatible_record_triples(profile)
        ),
        Fraction(0),
    )


def inserted_coefficients() -> dict[ProfileSplit, float]:
    cubic_target = occupation.endpoint_singleton_slice_energies(ORDER)[2] ** 0.5
    result = two_tier_coefficients(cubic_target, 0.5)
    for entry in block_coherent_entries():
        result[entry] = float(block_coherent_coefficient(*entry))
    return result


def diagnostic() -> BlockCoherentResult:
    entries = block_coherent_entries()
    exact = tuple(
        (block_coherent_coefficient(profile, split), profile, split)
        for profile, split in entries
    )
    minimum, _, _ = min(exact)
    maximum, maximum_profile, maximum_split = max(exact)
    record_sectors = sum(
        len(compatible_record_triples(profile)) for profile, _ in entries
    )
    cubic_target = occupation.endpoint_singleton_slice_energies(ORDER)[2] ** 0.5
    previous = optimize(mapped_coefficients=two_tier_coefficients(cubic_target, 0.5))
    inserted = optimize(mapped_coefficients=inserted_coefficients())
    return BlockCoherentResult(
        order=ORDER,
        dimension=DIMENSION,
        sign_modes=MODES,
        closed_entries=len(entries),
        record_sector_bounds=record_sectors,
        minimum_coefficient=float(minimum),
        minimum_coefficient_exact=str(minimum),
        maximum_coefficient=float(maximum),
        maximum_coefficient_exact=str(maximum),
        maximum_entry_profile=maximum_profile,
        maximum_entry_split=maximum_split,
        previous_two_tier=previous,
        block_coherent_inserted=inserted,
        margin_improvement=(
            inserted.margin_to_one_third - previous.margin_to_one_third
        ),
    )


def artifact_text(result: BlockCoherentResult) -> str:
    payload = {
        "schema": "round4_q64_block_coherent_contraction_v1",
        "result": asdict(result),
        "evidence_label": (
            "exact rational record-sector coefficients from the accepted "
            "weighted three-link path theorem; floating q64 Perron insertion; "
            "one batch only"
        ),
        "remaining_open": (
            "818 internally split balanced entries plus interval certification "
            "and the adaptive lift"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic q64 block-coherent result",
    )
    arguments = parser.parse_args()
    result = diagnostic()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "q64_block_coherent_contraction.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "q64 block-coherent contraction: "
        f"entries={result.closed_entries},"
        f"record_sectors={result.record_sector_bounds},"
        f"coefficient_min={result.minimum_coefficient:.12g},"
        f"coefficient_max={result.maximum_coefficient:.12g},"
        f"total={result.block_coherent_inserted.total:.12g},"
        f"margin={result.block_coherent_inserted.margin_to_one_third:.12g},"
        f"margin_gain={result.margin_improvement:.12g},"
        "status=proved_arbitrary_law_one_batch_entries"
    )


if __name__ == "__main__":
    main()
