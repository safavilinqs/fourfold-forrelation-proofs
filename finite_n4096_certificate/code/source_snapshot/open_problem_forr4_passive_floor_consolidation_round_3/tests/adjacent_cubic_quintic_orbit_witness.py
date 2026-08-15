#!/usr/bin/env python3
"""Regression for the adjacent cubic--quintic orbit and slice reductions."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from adjacent_cubic_quintic_orbit_witness import (
    direct_adjacent_orbit_witness,
    exact_link_moments,
    horizontal_adjacent_slice_certificate,
    horizontal_record_three_slice_sum,
    joint_record_one_slice_sum,
    parity_record_size,
    record_one_link_moment,
    record_three_link_moment,
    reduced_adjacent_orbit_witness,
)


def exact_q4_moment_checks() -> None:
    moments = exact_link_moments(4)
    cubic_record_one = [
        index
        for index, support in enumerate(moments.supports_three)
        if parity_record_size(4, support, axis=0) == 1
        and parity_record_size(4, support, axis=1) == 1
    ]
    quintic_record_one = [
        index
        for index, support in enumerate(moments.supports_five)
        if parity_record_size(4, support, axis=0) == 1
    ]
    maximum_error = 0.0
    for cubic_index in cubic_record_one:
        cubic = moments.supports_three[cubic_index]
        for quintic_index in quintic_record_one:
            maximum_error = max(
                maximum_error,
                abs(
                    record_one_link_moment(
                        4,
                        cubic,
                        moments.supports_five[quintic_index],
                    )
                    - moments.moment_35[cubic_index, quintic_index]
                ),
            )
    if maximum_error != 0:
        raise AssertionError(("record-one M35 formula", maximum_error))

    cubic_record_three = [
        index
        for index, support in enumerate(moments.supports_three)
        if parity_record_size(4, support, axis=0) == 1
        and parity_record_size(4, support, axis=1) == 3
    ]
    quintic_record_three = [
        index
        for index, support in enumerate(moments.supports_five)
        if parity_record_size(4, support, axis=0) == 3
    ]
    for cubic_index in cubic_record_three:
        cubic = moments.supports_three[cubic_index]
        for quintic_index in quintic_record_three:
            observed = record_three_link_moment(
                4,
                cubic,
                moments.supports_five[quintic_index],
            )
            expected = moments.moment_35[cubic_index, quintic_index]
            if observed != expected:
                raise AssertionError(
                    (
                        "record-three M35 formula",
                        cubic_index,
                        quintic_index,
                        observed,
                        expected,
                    )
                )


def main() -> None:
    exact_q4_moment_checks()
    direct = direct_adjacent_orbit_witness(4, 4, 1, (0, 4, 1))
    reduced = reduced_adjacent_orbit_witness(4, 4, 1, (0, 4, 1))
    if not np.isclose(direct.coefficient, 17 / 576, atol=2e-14):
        raise AssertionError(("direct q4 coefficient", direct))
    if not np.isclose(reduced.coefficient, direct.coefficient, atol=2e-14):
        raise AssertionError(("direct/reduced q4 coefficient", direct, reduced))
    if direct.rank != 16 * reduced.rank:
        raise AssertionError(("M11 rank multiplicity", direct.rank, reduced.rank))

    q4_formula = horizontal_adjacent_slice_certificate(4)
    q4_record_one = joint_record_one_slice_sum(4, 0, (0, 1, 2))
    q4_record_three = horizontal_record_three_slice_sum(4)
    if not np.isclose(
        q4_record_one.squared_moment_sum,
        q4_formula.record_one_squared_sum,
        atol=2e-13,
    ):
        raise AssertionError(("horizontal record-one formula", q4_record_one))
    if not np.isclose(
        q4_record_three.squared_moment_sum,
        q4_formula.record_three_squared_sum,
        atol=2e-13,
    ):
        raise AssertionError(("horizontal record-three formula", q4_record_three))

    q8 = reduced_adjacent_orbit_witness(8, 8, 1, (0, 16, 2))
    if not np.isclose(q8.coefficient, 2071 / 1204224, atol=2e-14):
        raise AssertionError(("q8 reduced orbit", q8))

    target = horizontal_adjacent_slice_certificate(32)
    if target.record_one_extensions != (1398, 43152):
        raise AssertionError(("q32 record-one counts", target))
    if not np.isclose(
        target.record_one_coefficient,
        0.0118822088518128,
        atol=2e-14,
    ):
        raise AssertionError(("q32 record-one coefficient", target))
    if not np.isclose(
        target.record_three_squared_sum,
        141,
        atol=1e-14,
    ):
        raise AssertionError(("q32 record-three slice", target))
    if not np.isclose(
        target.combined_coefficient,
        0.0168662459036314,
        atol=2e-14,
    ):
        raise AssertionError(("q32 combined slice", target))
    print(
        "adjacent cubic-quintic witness passed: "
        f"q4_orbit={direct.coefficient:.12g},"
        f"q8_orbit={q8.coefficient:.12g},"
        f"q32_record_one={target.record_one_coefficient:.12g},"
        f"q32_record_three={target.record_three_coefficient:.12g},"
        f"q32_combined={target.combined_coefficient:.12g}"
    )


if __name__ == "__main__":
    main()
