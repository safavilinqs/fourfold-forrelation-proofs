#!/usr/bin/env python3
"""Regression for the true terminal level-twelve sigma-one witness."""

from __future__ import annotations

from itertools import permutations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "searches"))

from terminal_interpolation_sigma_one_witness import (  # noqa: E402
    TRANSFERS,
    relaxed_star_boundary_degree,
    replay_terminal_witness,
    terminal_boundary_degree_cap,
)


EXPECTED_EDGES = {
    ((0, 0), (1, 0)),
    ((1, 3), (2, 2)),
    ((2, 1), (3, 1)),
    ((0, 1), (1, 3)),
    ((1, 0), (2, 0)),
    ((1, 2), (2, 1)),
    ((2, 0), (3, 0)),
    ((2, 2), (3, 2)),
    ((0, 3), (1, 2)),
}


def valid_order_count() -> int:
    """Replay every ordering and count the valid Stein partial orders."""

    target = replay_terminal_witness()
    valid = 0
    for order in permutations(TRANSFERS):
        try:
            witness = replay_terminal_witness(order)
        except ValueError:
            continue
        if witness.vertices != target.vertices or set(witness.edges) != EXPECTED_EDGES:
            raise AssertionError(("valid order changes the terminal graph", order))
        if not witness.local_weight_strictly_positive:
            raise AssertionError(("valid order loses positive local weight", order))
        valid += 1
    return valid


def main() -> None:
    witness = replay_terminal_witness()
    if len(witness.vertices) != 12 or len(witness.edges) != 9:
        raise AssertionError(("level-twelve profile", witness))
    if set(witness.edges) != EXPECTED_EDGES:
        raise AssertionError(("terminal edge set", witness.edges))
    if tuple(map(len, witness.components)) != (4, 4, 4):
        raise AssertionError(("three path components", witness.components))
    for component in witness.components:
        if {layer for layer, _ in component} != {0, 1, 2, 3}:
            raise AssertionError(("component misses a layer", component))
    if len(witness.edges) != len(witness.vertices) - len(witness.components):
        raise AssertionError(("terminal graph is not a forest", witness))
    if witness.initial_potential != 12 or witness.terminal_potential != 12:
        raise AssertionError(("all-new potential", witness))
    if witness.new_vertex_transfers != 6 or witness.existing_vertex_transfers != 0:
        raise AssertionError(("transfer classes", witness))
    if witness.first_layer_vertices != 3 or not witness.reflection_sensitive:
        raise AssertionError(("antisymmetrization sensitivity", witness))
    if witness.projective_sigma != 3 or witness.assigned_sigma != 1:
        raise AssertionError(("exact suppression parameters", witness))
    if not witness.local_weight_strictly_positive:
        raise AssertionError(("nonzero scalar branching coefficient", witness))
    if witness.maximum_boundary_degree > terminal_boundary_degree_cap():
        raise AssertionError(("terminal boundary-degree cap", witness))
    if not relaxed_star_boundary_degree(12) > terminal_boundary_degree_cap():
        raise AssertionError("old relaxed star should be outside the true image")

    # At level twelve, sigma=1 gives N^{-1/2} and therefore the unchanged
    # dose exponent 1/(2*12)=1/24.  The desired N^{1/16} row would require
    # integer sigma at least two.
    if 1 / (2 * 12) != 1 / 24:
        raise AssertionError("level-twelve exponent arithmetic")
    if witness.assigned_sigma >= 2:
        raise AssertionError("the structural exclusion route unexpectedly passed")

    order_count = valid_order_count()
    print(
        "terminal interpolation sigma-one witness passed: "
        f"v={len(witness.vertices)},"
        f"e={len(witness.edges)},"
        f"components={len(witness.components)},"
        f"valid_transfer_orders={order_count},"
        f"first_layer={witness.first_layer_vertices},"
        f"projective_sigma={witness.projective_sigma},"
        f"assigned_sigma={witness.assigned_sigma},"
        f"old_star_degree={relaxed_star_boundary_degree(12)},"
        f"true_boundary_cap={terminal_boundary_degree_cap()}"
    )


if __name__ == "__main__":
    main()
