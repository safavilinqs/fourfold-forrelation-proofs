#!/usr/bin/env python3
"""Protect the Round 4 mission, inherited baseline, and decision gates."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

REQUIRED_FILES = (
    "README.md",
    "GOAL.md",
    "JOINT_DECISION_TARGET.md",
    "INHERITED_BASELINE.md",
    "ROUND_4_CHARTER.md",
    "PLAN.md",
    "PAPER_TARGET.md",
    "EXPERIMENTAL_SCORECARD.md",
    "DECISION_GATES.md",
    "RESULTS_LEDGER.md",
    "OPEN_PROBLEMS.md",
    "notes/FIRST_PROJECT.md",
    "notes/Q64_DUAL_ENDPOINT_SCHUR_INSERTION.md",
    "notes/Q64_DUAL_ENDPOINT_INDEPENDENT_AUDIT.md",
    "notes/Q64_COMPLETE_OUTWARD_LEDGER.md",
    "notes/Q64_ADAPTIVE_TREE_FRONTIER_THEOREM.md",
    "notes/PAPER_CLOSEOUT_THEOREM_PACKAGE.md",
    "notes/EXPERIMENTAL_FEASIBILITY_DECISION.md",
    "notes/Q64_DECORATED_ADJACENT_ROW_INSERTION.md",
    "notes/Q64_DEGREE_TEN_COMPLETION_ROW_INSERTION.md",
    "notes/Q64_WHOLE_CUBIC_DECORATED_ROW_INSERTION.md",
    "notes/Q64_LAST_DEGREE_TEN_CHAIN_INSERTION.md",
    "notes/Q64_INTERNAL_WHOLE_CUBIC_ENDPOINT_INSERTION.md",
    "notes/Q64_BALANCED_PAIR_TRIPLE_MASK_INSERTION.md",
    "notes/Q64_ADJACENT_DOUBLE_CUBIC_QUINTIC_ENDPOINT_INSERTION.md",
    "notes/Q64_SHARED_QUINTIC_AND_ADAPTIVE_ACCEPTANCE.md",
    "notes/Q64_SHARED_QUINTIC_ROW_CHAIN_INSERTION.md",
    "notes/Q64_NONCUBIC_RECOVERED_UNIVERSAL_INSERTION.md",
    "notes/Q64_WHOLE_HIGHER_SPLIT_CUBIC_INSERTION.md",
    "notes/Q64_SAME_SIDE_WHOLE_LINK_INSERTION.md",
    "notes/Q64_MASKED_UNIVERSAL_AUDIT.md",
    "notes/Q64_MASKED_QUINTIC_SLICE_REPAIR.md",
    "notes/Q64_MASKED_LOCAL_WALSH_REPAIR.md",
    "notes/Q64_MASKED_CUBIC_ENDPOINT_REPAIR.md",
    "notes/Q64_MASKED_DOUBLE_QUINTIC_ENDPOINT_REPAIR.md",
    "notes/Q64_MASKED_DOUBLE_QUINTIC_RECORD_REPAIR.md",
    "notes/Q64_MASKED_FOUR_CUBIC_INCIDENCE_REPAIR.md",
    "notes/Q64_MASKED_CUBIC_SEPTIMIC_CHAIN_REPAIR.md",
    "notes/Q64_MASKED_RECOVERED_CUBIC_QUINTIC_INCIDENCE_REPAIR.md",
    "notes/Q64_RECOVERED_CUBIC_QUINTIC_INDEPENDENT_AUDIT.md",
    "notes/Q64_JOINT_RECOVERED_CUBIC_QUINTIC_CONTRACTION.md",
    "notes/Q64_FINAL_RESIDUAL_CHAIN_CONTRACTION.md",
    "notes/MASKED_TRANSLATION_REDUCTION.md",
    "artifacts/q64_dual_endpoint_schur_insertion.json",
    "artifacts/q64_masked_universal_audit.json",
    "artifacts/q64_masked_quintic_slice_repair.json",
    "artifacts/q64_masked_local_walsh_repair.json",
    "artifacts/q64_masked_cubic_endpoint_repair.json",
    "artifacts/q64_masked_double_quintic_endpoint_repair.json",
    "artifacts/q64_masked_double_quintic_record_repair.json",
    "artifacts/q64_masked_four_cubic_incidence_repair.json",
    "artifacts/q64_masked_cubic_septimic_chain_repair.json",
    "artifacts/q64_masked_recovered_cubic_quintic_incidence_repair.json",
    "artifacts/q64_joint_recovered_cubic_quintic_contraction.json",
    "artifacts/q64_final_residual_chain_contraction.json",
    "artifacts/q64_complete_outward_ledger.json",
    "artifacts/q64_adaptive_tree_frontier.json",
    "artifacts/q64_experimental_feasibility_gate.json",
    "artifacts/q8_masked_separated_quintic_residual_screen.json",
    "artifacts/q4_masked_translation_subspace_screen.json",
    "artifacts/q8_masked_translation_subspace_screen.json",
    "artifacts/q4_masked_translation_mixture_screen.json",
    "artifacts/q8_masked_translation_mixture_screen.json",
    "artifacts/q4_masked_translation_cocycle_inventory.json",
    "artifacts/q8_masked_translation_cocycle_inventory.json",
    "artifacts/q64_masked_translation_cocycle_inventory.json",
    "artifacts/q4_masked_translation_full_group_screen.json",
)


def require(path: str, *phrases: str) -> None:
    """Require every phrase in one Round 4 document."""

    text = (ROOT / path).read_text()
    for phrase in phrases:
        if phrase not in text:
            raise AssertionError((path, phrase))


def main() -> None:
    missing = [path for path in REQUIRED_FILES if not (ROOT / path).is_file()]
    if missing:
        raise AssertionError(("missing Round 4 files", missing))

    require(
        "GOAL.md",
        "distinctness-masked",
        "all 354",
        "888 balanced entries",
        "0.268858135059926",
        "passive hard dose",
        "greater than six",
        "active hard-dose-six protocol",
        "initial feasibility window",
        "experimentally credible",
        "branchwise hard photon-pass dose",
    )
    require(
        "JOINT_DECISION_TARGET.md",
        "historical initial decision framework",
        "leading 16 orbits",
        "leading 51 orbits",
        "90.5 percent",
        "proof architecture",
        "quadratic-bent",
    )
    require(
        "INHERITED_BASELINE.md",
        r"\Omega(N^{1/12})",
        "inherits no proof",
        "1,080 positive all-fresh histories",
        "848",
    )
    require(
        "ROUND_4_CHARTER.md",
        "Masked q64 contraction and registry",
        "Structural coefficient closure",
        "Dual caveat and ledger rebuild",
        "Adaptive viability and theorem",
        "Experimental resource statement",
    )
    require(
        "PLAN.md",
        "Phase 1",
        "masked q64 foundation",
        "compatible lower-witness test",
        "masked coefficient-one lemma",
        "smallest credible certified size",
        "experimental resource translation",
        "Explicitly deferred",
    )
    require(
        "PAPER_TARGET.md",
        r"D_{\mathsf P}^{\rm hard}(N_*)>6",
        r"N_*=4096",
        r"M=4N_*=16{,}384",
        "Kill condition",
        "not yet experimentally credible",
    )
    require(
        "EXPERIMENTAL_SCORECARD.md",
        "active sample traversals",
        "passive dose excluded",
        "loss model",
        "detector assumptions",
        "VERDICT: NOT YET EXPERIMENTALLY CREDIBLE",
    )
    require(
        "DECISION_GATES.md",
        "Diagnostics and relaxations must never be added",
        "masked coefficient-one decision",
        "all 80 later entries",
        "remaining 12 caveated entries",
        "0.063475198273408",
        "experimental credibility",
        "paper package",
        "Do not reopen terminal-shape enumeration",
    )
    require(
        "RESULTS_LEDGER.md",
        "Paper headline target",
        "Passive lower at",
        "Active upper",
        "experimental feasibility",
        "unresolved complement/reversal orbits",
        "Foundational masked-universal audit",
        "Historical conditional noncubic",
        "Preserved same-side whole-link theorem",
    )
    require(
        "notes/FIRST_PROJECT.md",
        "historical and superseded",
        "Deliverable 1",
        "compatible lower-witness vector",
        "shared upper contraction",
        "Pass",
        "Fail",
        "Deliverable 6",
        "Stop rules",
    )
    require(
        "notes/Q64_DUAL_ENDPOINT_SCHUR_INSERTION.md",
        "twelve entries",
        "0.149556115743",
        "Schur-multiplier composition",
        "888 of 888",
    )
    require(
        "notes/Q64_DECORATED_ADJACENT_ROW_INSERTION.md",
        "16 entries",
        "0.0200795672469",
        "unsplit outer cubic",
        "136 entries",
    )
    require(
        "notes/Q64_DEGREE_TEN_COMPLETION_ROW_INSERTION.md",
        "12 entries",
        "0.00861554231015",
        "0.0754041939294",
        "124 entries",
    )
    require(
        "notes/Q64_WHOLE_CUBIC_DECORATED_ROW_INSERTION.md",
        "16 entries",
        "0.000196155632204",
        "0.213518291069",
        "108 entries",
    )
    require(
        "notes/Q64_LAST_DEGREE_TEN_CHAIN_INSERTION.md",
        "final four-entry degree-ten orbit",
        "0.0910312181521",
        "104 such entries",
        "1/q",
    )
    require(
        "notes/Q64_INTERNAL_WHOLE_CUBIC_ENDPOINT_INSERTION.md",
        "16 degree-twelve entries",
        "0.113036239514",
        "0.281119075921",
        "remaining quintic inventory is 88 entries",
    )
    require(
        "notes/Q64_BALANCED_PAIR_TRIPLE_MASK_INSERTION.md",
        "all eight",
        "1+\\sqrt6+\\sqrt3",
        "0.642693497508",
        "80 quintic entries",
    )
    require(
        "notes/Q64_ADJACENT_DOUBLE_CUBIC_QUINTIC_ENDPOINT_INSERTION.md",
        "eight four-entry",
        "4.70180143564",
        "48 quintic entries",
        "incidence theorem was proved",
    )
    require(
        "notes/Q64_SHARED_QUINTIC_AND_ADAPTIVE_ACCEPTANCE.md",
        "48 entries",
        "0.410314553367",
        "0.008521770161998",
        "outcome width",
        "Kill and pivot conditions",
    )
    require(
        "notes/Q64_SHARED_QUINTIC_ROW_CHAIN_INSERTION.md",
        "all 48 previously open",
        "0.0203737451368",
        "five finite row/chain templates",
        "0.008521770161998",
        "remaining 460 entries",
    )
    require(
        "notes/Q64_NONCUBIC_RECOVERED_UNIVERSAL_INSERTION.md",
        "mixed and partially quarantined",
        "0.0382305883153",
        "96-entry two-split-cubic",
        "664-entry count",
        "values are withdrawn",
    )
    require(
        "notes/Q64_WHOLE_HIGHER_SPLIT_CUBIC_INSERTION.md",
        "48-entry mask-aware theorem is preserved",
        "0.00846466875312",
        "cumulative 712-entry count",
        "only 176 entries remain",
        "0.142581909211",
    )
    require(
        "notes/Q64_SAME_SIDE_WHOLE_LINK_INSERTION.md",
        "96-entry mask-aware theorem is preserved",
        "cumulative 808-entry count",
        "0.167292848473",
    )
    require(
        "notes/Q64_MASKED_UNIVERSAL_AUDIT.md",
        "180 local-Walsh",
        "12 cubic-endpoint",
        "6 double-quintic endpoint",
        "38 four-cubic incidence",
        "12 cubic--septimic chain",
        "28 recovered cubic--quintic endpoint-row",
        "joint shared-quintic",
        "12 double-quintic record",
        "888 certified",
        "zero open",
        "complete physical signed-permutation occurrence matrix",
        "3.40344112205",
    )
    require(
        "notes/Q64_MASKED_QUINTIC_SLICE_REPAIR.md",
        "proved arbitrary-correlated-diagonal",
        "54 entries",
        "actual distinctness-masked occurrence matrix",
        "300",
        "0.646562122163",
    )
    require(
        "notes/Q64_MASKED_LOCAL_WALSH_REPAIR.md",
        "theorem for 180 actual",
        "120 of the original 354",
        "54\\over64",
        "889856 exact comparisons",
    )
    require(
        "notes/Q64_MASKED_CUBIC_ENDPOINT_REPAIR.md",
        "theorem for 12 actual",
        "108 remain quarantined",
        "225\\over4",
        "1920 physical rows",
    )
    require(
        "notes/Q64_MASKED_DOUBLE_QUINTIC_ENDPOINT_REPAIR.md",
        "theorem for six actual masked",
        "252 of 354 affected",
        "102 remain quarantined",
        "Schur-multiplier factorization",
        "0.999511599483",
    )
    require(
        "notes/Q64_MASKED_DOUBLE_QUINTIC_RECORD_REPAIR.md",
        "all twelve residual entries",
        "342 entries",
        "12 recovered",
        "5901977909483",
        "0.00213806918675",
        "At $q=8$",
    )
    require(
        "notes/Q64_MASKED_FOUR_CUBIC_INCIDENCE_REPAIR.md",
        "theorem for all 38 residual",
        "290 of the 354",
        "64 remain quarantined",
        "157,952 compatible",
        "0.00894228260682",
    )
    require(
        "notes/Q64_MASKED_CUBIC_SEPTIMIC_CHAIN_REPAIR.md",
        "theorem for all twelve",
        "302",
        "52 remain quarantined",
        "110311919",
        "0.0884219637995",
        "0.538626546351",
    )
    require(
        "notes/Q64_MASKED_RECOVERED_CUBIC_QUINTIC_INCIDENCE_REPAIR.md",
        "theorem for 28 recovered",
        "twelve $(1,3,5,3)$",
        "1/17920",
        "0.0162888571820",
        "0.703615181088",
    )
    require(
        "notes/Q64_RECOVERED_CUBIC_QUINTIC_INDEPENDENT_AUDIT.md",
        "VERDICT: REJECTED",
        "1/17920",
        "1/25088",
        "342 of 354",
        "joint endpoint-compatible contraction",
    )
    require(
        "notes/Q64_JOINT_RECOVERED_CUBIC_QUINTIC_CONTRACTION.md",
        "VERDICT: CERTIFIED",
        "15 feasible joint shapes",
        "0.33828697324447987",
        r"\frac1{17920}",
        r"\frac1{25088}",
        "354 of 354",
    )
    require(
        "notes/Q64_FINAL_RESIDUAL_CHAIN_CONTRACTION.md",
        "VERDICT: CERTIFIED",
        "All eighty entries",
        "forty-eight",
        "thirty-two",
        "18816",
        "0.447299757774",
        "arbitrary correlated",
    )
    require(
        "notes/Q64_DUAL_ENDPOINT_INDEPENDENT_AUDIT.md",
        "VERDICT: CERTIFIED",
        r"\frac{1985}{129024}",
        r"\frac{62527}{43008}",
        "192 exact completed-link Gram identities",
        "888 of 888",
    )
    require(
        "notes/Q64_COMPLETE_OUTWARD_LEDGER.md",
        "6{,}016",
        "888",
        "0.268858135059925",
        "0.063475198273408",
        "Collatz--Wielandt",
        "tree-frontier theorem",
    )
    require(
        "notes/Q64_ADAPTIVE_TREE_FRONTIER_THEOREM.md",
        "VERDICT: CERTIFIED",
        "multiplier is exactly one",
        "0.268858135059925",
        "0.063475198273408",
        "direct-sum features",
        "Rare outcomes",
        "classical feed-forward",
    )
    require(
        "notes/MASKED_TRANSLATION_REDUCTION.md",
        "97 disjoint structural orbits",
        "arbitrary-law twirling reduction",
        "0.0220970869121",
        "0.000279017854455",
        "mixture of orbit shapes",
        "Exact projective-type classification",
        "| 0 | 21 |",
        "| 4 | 26 |",
        "| 8 | 50 |",
        "Exact full-group Clifford formula",
        "0.176776695297",
        "shape-indexed symbol",
    )
    require(
        "notes/FINITE_SIZE_LEDGER_AUDIT.md",
        "888 entries",
        "224 complement/reversal orbits",
        "844 entries",
        "0.041462318296515",
        r"N\in\{256,1024,4096\}",
    )

    print(
        "round-four initialization passed: "
        f"files={len(REQUIRED_FILES)},"
        "primary=q64_masked_coefficient_one_decision,"
        "lead_N=4096,"
        "certified_entries=888,"
        "conservative_supported_entries=888,"
        "quarantined_entries=0,"
        "conservative_open_entries=0,"
        "outward_ledger=passes_reserve,"
        "adaptive_frontier_multiplier=1,"
        "active_dose=6,"
        "asymptotic_baseline=1/12,"
        "sharp_grouped_barrier=preserved"
    )


if __name__ == "__main__":
    main()
