#!/usr/bin/env python3
"""Machine-readable active six-dose resource row at the q64 target size."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from json import dumps
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DIMENSION = 4096
SIGN_BLOCKS = 4
SIGN_MODES = SIGN_BLOCKS * DIMENSION
FLAGS = 3
HARD_DOSE_PER_FLAG = 2


@dataclass(frozen=True)
class BranchWord:
    path: str
    chronological_operations: tuple[str, ...]
    charged_sign_blocks: tuple[int, ...]
    public_hadamards: int
    hard_dose: int


@dataclass(frozen=True)
class ActiveResourceRow:
    dimension: int
    sign_blocks: int
    sign_modes: int
    flags: int
    single_photons: int
    logical_path_dimension_per_flag: int
    logical_mode_dimension_per_flag: int
    hard_dose_per_flag: int
    total_hard_dose: int
    left_branch: BranchWord
    right_branch: BranchWord
    state_preparation: str
    receiver: str
    postselection: bool
    flag_correct_probability_at_promise_boundary: str
    majority_error_exact: str
    majority_error: float
    margin_below_one_third_exact: str
    margin_below_one_third: float


def resource_row() -> ActiveResourceRow:
    flag_correct = Fraction(5, 8)
    flag_error = 1 - flag_correct
    majority_error = flag_error**3 + 3 * flag_correct * flag_error**2
    margin = Fraction(1, 3) - majority_error
    left = BranchWord(
        path="left",
        chronological_operations=("D1", "H", "D2"),
        charged_sign_blocks=(1, 2),
        public_hadamards=1,
        hard_dose=2,
    )
    right = BranchWord(
        path="right",
        chronological_operations=("D4", "H", "D3", "H"),
        charged_sign_blocks=(4, 3),
        public_hadamards=2,
        hard_dose=2,
    )
    return ActiveResourceRow(
        dimension=DIMENSION,
        sign_blocks=SIGN_BLOCKS,
        sign_modes=SIGN_MODES,
        flags=FLAGS,
        single_photons=FLAGS,
        logical_path_dimension_per_flag=2,
        logical_mode_dimension_per_flag=DIMENSION,
        hard_dose_per_flag=HARD_DOSE_PER_FLAG,
        total_hard_dose=FLAGS * HARD_DOSE_PER_FLAG,
        left_branch=left,
        right_branch=right,
        state_preparation=(
            "one photon in |+>_path tensor |u>_mode for each independent flag"
        ),
        receiver=(
            "Pauli-X path interferometer with mode-insensitive port detection; "
            "majority vote over three retained binary flags"
        ),
        postselection=False,
        flag_correct_probability_at_promise_boundary=str(flag_correct),
        majority_error_exact=str(majority_error),
        majority_error=float(majority_error),
        margin_below_one_third_exact=str(margin),
        margin_below_one_third=float(margin),
    )


def artifact_text(result: ActiveResourceRow) -> str:
    payload = {
        "schema": "round4_active_six_resource_row_v1",
        "result": asdict(result),
        "theorem_scope": (
            "pointwise for every promised input at every Sylvester dimension; "
            "deterministic branchwise hard dose; no postselection"
        ),
        "experimental_scope": (
            "ideal coherent signs, public Hadamards, routing, source, and "
            "detection; loss and phase robustness not yet quantified"
        ),
    }
    return dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help="write the deterministic active resource row under artifacts/",
    )
    arguments = parser.parse_args()
    result = resource_row()
    if arguments.write_artifact:
        path = ROOT / "artifacts" / "active_six_resource_row.json"
        path.write_text(artifact_text(result), encoding="utf-8")
    print(
        "active six-dose resource row: "
        f"N={result.dimension},"
        f"M={result.sign_modes},"
        f"flags={result.flags},"
        f"photons={result.single_photons},"
        f"dose_per_flag={result.hard_dose_per_flag},"
        f"total_dose={result.total_hard_dose},"
        f"error={result.majority_error_exact},"
        f"margin={result.margin_below_one_third_exact},"
        "status=proved_ideal_active_protocol"
    )


if __name__ == "__main__":
    main()
