#!/usr/bin/env python3
"""Rigorous coarse frontier for the open degree-ten/twelve sectors.

The remaining high-degree support profiles are split by the three odd-label
record sizes carried by the links.  Two elementary bounds are then combined:

* a record sector has maximum entry at most
  ``prod(1 / comb(q, record))``;
* rank--Frobenius is sharpened by the maximum number of compatible support
  extensions of either a fixed row support or a fixed column support.

The one-axis extension count is exact.  For a block constrained on both axes,
the minimum of the two one-axis maxima is a safe relaxation.  The script also
applies one rank--incidence estimate directly to the whole profile, using the
universal maximum-entry bound ``q^-3``.  Taking the better of this estimate
and the triangle sum over record sectors gives a proved coefficient for every
open profile/split.  It is deliberately coarse: entries above ``1/q`` remain
open rather than being treated as counterexamples.

Finally, the exact occupation pairing filters the inventory to the balanced
splits that can actually join two total-occupation-six states.  A Perron
sensitivity calculation ranks those remaining symmetry orbits without
promoting the diagnostic target coefficients to theorems.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from math import comb, prod, sqrt

from occupation_compatible_sector_optimization import (
    HIGH_DEGREE_PROFILES,
    KNOWN_HIGH_DEGREE_PROFILES,
    certificate,
    integer_partitions,
    multiplicity,
    occupation_states,
    paired_state,
    profile_splits,
)


ORDER = 32
TARGET = 1 / ORDER

Profile = tuple[int, ...]
Split = tuple[int, ...]
Records = tuple[int, int, int]
ProfileSplit = tuple[Profile, Split]


@dataclass(frozen=True)
class FrontierSummary:
    profiles: int
    profile_splits: int
    record_sectors: int
    certified_profile_splits: int
    unresolved_profile_splits: int
    certified_degree_ten: int
    degree_ten_splits: int
    certified_degree_twelve: int
    degree_twelve_splits: int
    dose_six_relevant_splits: int
    certified_dose_six_relevant_splits: int


@dataclass(frozen=True)
class RoutePriority:
    entries: tuple[ProfileSplit, ...]
    perron_sensitivity: float
    perron_contribution: float
    current_coefficients: tuple[float, ...]


def open_profiles() -> tuple[Profile, ...]:
    """Return the degree-ten/twelve profiles not in the accepted ledger."""

    return tuple(
        profile
        for profile in HIGH_DEGREE_PROFILES
        if profile not in KNOWN_HIGH_DEGREE_PROFILES
    )


def compatible_record_triples(profile: Profile) -> tuple[Records, ...]:
    """Return all parity-compatible odd record sizes on the three links."""

    choices = tuple(
        range(1, min(left, right) + 1, 2)
        for left, right in zip(profile[:-1], profile[1:], strict=True)
    )
    return tuple(product(*choices))


def record_extension_count(
    order: int,
    fixed_axis_counts: tuple[int, ...],
    degree: int,
    record: int,
) -> int:
    """Count extensions with a prescribed odd-label record on one axis.

    ``fixed_axis_counts`` gives the multiplicities of a fixed selected support
    on the occupied axis labels.  Additional cells are chosen independently
    within labels, while a two-coordinate dynamic program tracks total degree
    and the number of labels with odd final multiplicity.
    """

    fixed_size = sum(fixed_axis_counts)
    if fixed_size > degree:
        return 0
    counts = fixed_axis_counts + (0,) * (order - len(fixed_axis_counts))
    dynamic = {(0, 0): 1}
    additions = degree - fixed_size
    for fixed in counts:
        updated: dict[tuple[int, int], int] = {}
        for (added, odd), count in dynamic.items():
            remaining = additions - added
            for extra in range(min(order - fixed, remaining) + 1):
                key = (added + extra, odd + (fixed + extra) % 2)
                updated[key] = updated.get(key, 0) + count * comb(
                    order - fixed, extra
                )
        dynamic = updated
    return dynamic.get((additions, record), 0)


@lru_cache(maxsize=None)
def record_incidence_for_order(
    order: int,
    degree: int,
    selected: int,
    record: int,
) -> int:
    """Maximum one-axis record-sector incidence over selected supports."""

    return max(
        record_extension_count(order, partition, degree, record)
        for partition in integer_partitions(selected)
        if len(partition) <= order and max(partition, default=0) <= order
    )


def block_record_incidence(
    order: int,
    degree: int,
    selected: int,
    left_record: int | None,
    right_record: int | None,
) -> int:
    """Safe incidence for a block with one or two record constraints."""

    available = tuple(
        record
        for record in (left_record, right_record)
        if record is not None
    )
    return min(
        record_incidence_for_order(order, degree, selected, record)
        for record in available
    )


def record_entry_bound(order: int, records: Records) -> float:
    """Universal maximum moment entry in one record sector."""

    return 1 / prod(comb(order, record) for record in records)


def rank_factor(order: int, profile: Profile, split: Split) -> int:
    """Square root of the smaller crude matrix dimension."""

    selected = sum(split)
    smaller_side = min(selected, sum(profile) - selected)
    return order**smaller_side


def sector_coefficient_bound(
    order: int,
    profile: Profile,
    split: Split,
    records: Records,
) -> float:
    """Rank--incidence Frobenius bound for one record sector."""

    left_records = (None,) + records
    right_records = records + (None,)
    row_incidence = prod(
        block_record_incidence(
            order,
            degree,
            selected,
            left,
            right,
        )
        for degree, selected, left, right in zip(
            profile,
            split,
            left_records,
            right_records,
            strict=True,
        )
    )
    column_incidence = prod(
        block_record_incidence(
            order,
            degree,
            degree - selected,
            left,
            right,
        )
        for degree, selected, left, right in zip(
            profile,
            split,
            left_records,
            right_records,
            strict=True,
        )
    )
    geometry = min(
        rank_factor(order, profile, split),
        sqrt(row_incidence),
        sqrt(column_incidence),
    )
    return record_entry_bound(order, records) * geometry


def all_support_incidence(
    order: int, degree: int, selected: int
) -> int:
    """Number of degree-sized supports extending a fixed selected support."""

    return comb(order**2 - selected, degree - selected)


def whole_profile_coefficient_bound(
    order: int, profile: Profile, split: Split
) -> float:
    """Single rank--incidence bound before splitting by record sizes."""

    row_incidence = prod(
        all_support_incidence(order, degree, selected)
        for degree, selected in zip(profile, split, strict=True)
    )
    column_incidence = prod(
        all_support_incidence(order, degree, degree - selected)
        for degree, selected in zip(profile, split, strict=True)
    )
    geometry = min(
        rank_factor(order, profile, split),
        sqrt(row_incidence),
        sqrt(column_incidence),
    )
    return geometry / order**3


def triangle_sector_coefficient_bound(
    order: int, profile: Profile, split: Split
) -> float:
    """Safe triangle sum of the separate record-sector estimates."""

    return sum(
        sector_coefficient_bound(order, profile, split, records)
        for records in compatible_record_triples(profile)
    )


def combined_coefficient_bound(
    order: int, profile: Profile, split: Split
) -> float:
    """Best of the whole-profile and record-sector triangle estimates."""

    return min(
        whole_profile_coefficient_bound(order, profile, split),
        triangle_sector_coefficient_bound(order, profile, split),
    )


def certified_coefficients(
    order: int = ORDER,
) -> dict[ProfileSplit, float]:
    """Return open splits whose proved coarse coefficient is at most 1/q."""

    target = 1 / order
    result: dict[ProfileSplit, float] = {}
    for profile in open_profiles():
        for split in profile_splits(profile):
            coefficient = combined_coefficient_bound(order, profile, split)
            if coefficient <= target:
                result[(profile, split)] = coefficient
    return result


def dose_six_relevant_entries() -> tuple[ProfileSplit, ...]:
    """Return splits that can join two total-occupation-six states.

    If the profile degree is ``L`` and the selected side has size ``k``,
    occupation compatibility sends total occupation six to ``6 + L - 2k``.
    Both sides can therefore have total six only when ``k=L/2``.
    """

    return tuple(
        (profile, split)
        for profile in open_profiles()
        for split in profile_splits(profile)
        if 2 * sum(split) == sum(profile)
    )


def unresolved_entries(order: int = ORDER) -> tuple[ProfileSplit, ...]:
    """Return entries not yet certified at the common 1/q target."""

    certified = certified_coefficients(order)
    return tuple(
        (profile, split)
        for profile in open_profiles()
        for split in profile_splits(profile)
        if (profile, split) not in certified
    )


def frontier_summary(order: int = ORDER) -> FrontierSummary:
    """Summarize the proved and unresolved portions of the frontier."""

    profiles = open_profiles()
    certified = certified_coefficients(order)
    degree_totals: dict[int, int] = {10: 0, 12: 0}
    degree_certified: dict[int, int] = {10: 0, 12: 0}
    record_sectors = 0
    for profile in profiles:
        splits = profile_splits(profile)
        degree = sum(profile)
        degree_totals[degree] += len(splits)
        degree_certified[degree] += sum(
            (profile, split) in certified for split in splits
        )
        record_sectors += len(splits) * len(
            compatible_record_triples(profile)
        )
    total_splits = sum(degree_totals.values())
    relevant = dose_six_relevant_entries()
    return FrontierSummary(
        profiles=len(profiles),
        profile_splits=total_splits,
        record_sectors=record_sectors,
        certified_profile_splits=len(certified),
        unresolved_profile_splits=total_splits - len(certified),
        certified_degree_ten=degree_certified[10],
        degree_ten_splits=degree_totals[10],
        certified_degree_twelve=degree_certified[12],
        degree_twelve_splits=degree_totals[12],
        dose_six_relevant_splits=len(relevant),
        certified_dose_six_relevant_splits=sum(
            entry in certified for entry in relevant
        ),
    )


def leading_unresolved(
    order: int = ORDER, limit: int = 20
) -> tuple[tuple[float, Profile, Split], ...]:
    """Rank unresolved entries by their current coarse upper coefficient."""

    entries = (
        (combined_coefficient_bound(order, profile, split), profile, split)
        for profile, split in unresolved_entries(order)
    )
    return tuple(sorted(entries, reverse=True)[:limit])


def symmetry_orbits(
    entries: tuple[ProfileSplit, ...],
) -> tuple[tuple[ProfileSplit, ...], ...]:
    """Group entries under support complementation and path reversal."""

    available = set(entries)
    seen: set[ProfileSplit] = set()
    result = []
    for entry in entries:
        if entry in seen:
            continue
        profile, split = entry
        complement = tuple(
            degree - selected
            for degree, selected in zip(profile, split, strict=True)
        )
        reversed_profile = tuple(reversed(profile))
        reversed_split = tuple(reversed(split))
        reversed_complement = tuple(reversed(complement))
        orbit = tuple(
            sorted(
                {
                    entry,
                    (profile, complement),
                    (reversed_profile, reversed_split),
                    (reversed_profile, reversed_complement),
                }
                & available
            )
        )
        seen.update(orbit)
        result.append(orbit)
    return tuple(result)


def split_perron_sensitivity(
    profile: Profile,
    split: Split,
    beta: float,
    weights: dict[tuple[int, ...], float],
) -> float:
    """Derivative of the Perron objective with respect to one coefficient."""

    states = occupation_states()
    state_set = set(states)
    complement = tuple(
        degree - selected
        for degree, selected in zip(profile, split, strict=True)
    )
    result = 0.0
    for state in states:
        if any(
            occupation < selected
            for occupation, selected in zip(state, split, strict=True)
        ):
            continue
        partner = paired_state(state, profile, split)
        if partner not in state_set:
            continue
        left_count = multiplicity(state, split)
        right_count = multiplicity(partner, complement)
        result += sqrt(
            left_count
            * right_count
            * weights[state]
            * weights[partner]
        )
    return beta ** sum(profile) * result


def coarse_target_priorities(
    limit: int | None = 20,
) -> tuple[RoutePriority, ...]:
    """Rank relevant symmetry orbits in the current all-open target ledger.

    The ranking is a route-selection diagnostic, not a coefficient bound.  It
    decomposes the Perron objective at the optimized all-open target into the
    contributions of the still-unproved balanced profile/split orbits.
    """

    from repaired_open_profile_budget import (
        coarse_open_completion_coefficients,
        coarse_open_completion_target,
    )

    target = coarse_open_completion_target()
    coefficients = coarse_open_completion_coefficients()
    ledger = certificate(
        beta=target.optimal_beta,
        profile_split_coefficients=coefficients,
    )
    weights = dict(ledger.occupation_weights)
    sensitivities = {
        entry: split_perron_sensitivity(
            entry[0], entry[1], target.optimal_beta, weights
        )
        for entry in dose_six_relevant_entries()
    }
    priorities = tuple(
        sorted(
            (
                RoutePriority(
                    entries=orbit,
                    perron_sensitivity=sum(
                        sensitivities[entry] for entry in orbit
                    ),
                    perron_contribution=sum(
                        sensitivities[entry] * coefficients[entry]
                        for entry in orbit
                    ),
                    current_coefficients=tuple(
                        coefficients[entry] for entry in orbit
                    ),
                )
                for orbit in symmetry_orbits(dose_six_relevant_entries())
            ),
            key=lambda item: item.perron_contribution,
            reverse=True,
        )
    )
    return priorities if limit is None else priorities[:limit]


def main() -> None:
    summary = frontier_summary()
    print(f"order={ORDER}, target={TARGET:.12g}")
    print(
        "inventory: "
        f"profiles={summary.profiles}, "
        f"profile_splits={summary.profile_splits}, "
        f"record_sectors={summary.record_sectors}"
    )
    print(
        "certified at 1/q: "
        f"total={summary.certified_profile_splits}/"
        f"{summary.profile_splits}, "
        f"degree10={summary.certified_degree_ten}/"
        f"{summary.degree_ten_splits}, "
        f"degree12={summary.certified_degree_twelve}/"
        f"{summary.degree_twelve_splits}"
    )
    print(
        "dose-six relevant: "
        f"certified={summary.certified_dose_six_relevant_splits}/"
        f"{summary.dose_six_relevant_splits}"
    )
    print(f"unresolved={summary.unresolved_profile_splits}")
    for priority in coarse_target_priorities(limit=12):
        print(
            "dose-six route priority: "
            f"contribution={priority.perron_contribution:.12g}, "
            f"sensitivity={priority.perron_sensitivity:.12g}, "
            f"entries={priority.entries}"
        )


if __name__ == "__main__":
    main()
