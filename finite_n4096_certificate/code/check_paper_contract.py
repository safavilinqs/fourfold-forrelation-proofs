#!/usr/bin/env python3
"""Check manuscript constants, citation resolution, figures, and package structure."""

from __future__ import annotations

import json
import re
from decimal import Decimal
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
ROUND4 = (
    PROJECT
    / "code"
    / "source_snapshot"
    / "open_problem_forr4_passive_floor_consolidation_round_4"
)
ARTIFACTS = ROUND4 / "artifacts"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def load_json(name: str) -> dict:
    return json.loads((ARTIFACTS / name).read_text(encoding="utf-8"))


def manuscript_text() -> str:
    paths = [PROJECT / "main.tex", PROJECT / "macros.tex"]
    paths.extend(sorted((PROJECT / "sections").glob("*.tex")))
    paths.extend(sorted((PROJECT / "appendix").glob("*.tex")))
    return "\n".join(path.read_text(encoding="utf-8") for path in paths)


def check_structure() -> None:
    required = (
        "main.tex",
        "macros.tex",
        "references.bib",
        "README.md",
        "GOAL.md",
        "STATUS.md",
        "RUN_OF_RECORD.md",
        "dev/BIBLIOGRAPHY_NOTES.md",
        "dev/PAPER_PLAN.md",
        "dev/SOURCE_MAP.md",
        "dev/STYLE_CONTRACT.md",
        "figures/generate_figures.py",
        "figures/active_protocol.pdf",
        "figures/signed_permutation_phase_grids.pdf",
        "figures/passive_certificate.pdf",
    )
    for relative in required:
        require((PROJECT / relative).is_file(), f"missing required file: {relative}")
    for round_name in (2, 3, 4):
        folder = ROUND4.parent / (
            "open_problem_forr4_passive_floor_consolidation_round_"
            f"{round_name}"
        )
        require(folder.is_dir(), f"missing source snapshot: {folder.name}")


def check_artifacts(text: str) -> None:
    ledger = load_json("q64_complete_outward_ledger.json")["result"]
    adaptive = load_json("q64_adaptive_tree_frontier.json")["result"]
    active = load_json("active_six_resource_row.json")["result"]
    robustness = load_json("active_six_robustness_gate.json")["result"]
    feasibility = load_json("q64_experimental_feasibility_gate.json")

    require(ledger["dimension"] == 4096, "ledger dimension changed")
    require(
        ledger["certified_balanced_high_sector_coefficients"] == 888,
        "balanced theorem coefficient count changed",
    )
    require(ledger["supported_balanced_entries"] == 888, "registry count changed")
    require(ledger["open_balanced_entries"] == 0, "registry has open entries")
    require(
        ledger["excluded_unbalanced_high_sector_entries"] == 5128,
        "unbalanced routing inventory changed",
    )
    require(
        ledger["excluded_unbalanced_high_sector_incidence_records"] == 272,
        "unbalanced incidence audit changed",
    )
    require(
        ledger["excluded_unbalanced_high_sector_undirected_edges"] == 136,
        "unbalanced edge audit changed",
    )
    require(ledger["passes_reserve_gate"] is True, "ledger reserve gate failed")
    require(adaptive["frontier_multiplier_numerator"] == 1, "adaptive numerator changed")
    require(adaptive["frontier_multiplier_denominator"] == 1, "adaptive denominator changed")
    require(active["total_hard_dose"] == 6, "active dose changed")
    require(active["majority_error_exact"] == "81/256", "active error changed")
    require(feasibility["verdict"] == "NOT_YET_EXPERIMENTALLY_CREDIBLE", "feasibility verdict changed")

    displayed_uppers = {
        "total_upper": Decimal("0.260969224792207925"),
        "collatz_perron_upper": Decimal("0.258744096385577223"),
        "promise_loss_upper": Decimal("0.002225128406630703"),
    }
    for key, value in displayed_uppers.items():
        require(str(value) in text, f"manuscript does not contain displayed upper: {value}")
        require(value >= Decimal(str(ledger[key])), f"displayed {key} is rounded inward")

    displayed_reserve_lower = Decimal("0.071364108541125408")
    displayed_error_lower = Decimal("0.369515387603896037")
    exact_error_lower = (Decimal(1) - Decimal(str(ledger["total_upper"]))) / 2
    require(str(displayed_reserve_lower) in text, "reserve lower missing")
    require(str(displayed_error_lower) in text, "passive error lower missing")
    require(displayed_reserve_lower <= Decimal(str(ledger["reserve_margin_lower"])), "reserve lower rounded upward")
    require(displayed_error_lower <= exact_error_lower, "passive error lower rounded upward")

    expected_strings = (
        "888/888",
        "81/256",
        "0.904294855157",
        "total dose at most six",
        "block diagonal in total signal photon number",
        "272 unbalanced",
        "136 undirected edges",
        "not covered",
    )
    for value in expected_strings:
        require(value in text, f"manuscript does not contain required value: {value}")
    require("N^{1/8}" in text, "corrected asymptotic exponent is missing")
    require("N^{1/12}" not in text, "superseded asymptotic exponent remains")

    displayed_contrast = float("0.904294855157")
    require(
        displayed_contrast >= robustness["minimum_multiplicative_contrast"],
        "displayed contrast threshold is rounded inward",
    )


def check_citations(text: str) -> None:
    require("\\citetag" not in text, "literal AuthorYear citation macro found")
    require("\\nocite" not in text, "forced bibliography inclusion found")
    citation_groups = re.findall(r"\\cite(?:p|t)?\{([^}]+)\}", text)
    cited_keys = {
        key.strip()
        for group in citation_groups
        for key in group.split(",")
        if key.strip()
    }
    notes = (PROJECT / "dev" / "BIBLIOGRAPHY_NOTES.md").read_text(encoding="utf-8")
    bib = (PROJECT / "references.bib").read_text(encoding="utf-8")
    bib_keys = set(re.findall(r"@\w+\{([^,]+),", bib))
    require(cited_keys, "no conventional bibliography citations found")
    for key in sorted(cited_keys):
        require(f"## `{key}`" in notes, f"citation note missing: {key}")
        require(key in bib_keys, f"BibTeX entry missing: {key}")
    require(
        bib_keys == cited_keys,
        "bibliography contains uncited or cited-but-missing entries: "
        f"uncited={sorted(bib_keys - cited_keys)}, missing={sorted(cited_keys - bib_keys)}",
    )


def check_figures(text: str) -> None:
    stems = (
        "active_protocol",
        "signed_permutation_phase_grids",
        "passive_certificate",
    )
    for stem in stems:
        for suffix in ("pdf", "svg"):
            path = PROJECT / "figures" / f"{stem}.{suffix}"
            require(path.is_file() and path.stat().st_size > 0, f"missing figure: {path.name}")
        require(
            f"figures/{stem}.pdf" in text,
            f"figure is not included in manuscript: {stem}",
        )


def check_style(text: str) -> None:
    banned = (
        "It's not ",
        "It is not just ",
        "game-changer",
        "paradigm shift",
        "low-hanging fruit",
        "north star",
    )
    for phrase in banned:
        require(phrase not in text, f"banned prose pattern found: {phrase}")
    require("`" not in text, "Markdown backtick found in LaTeX source; use \\path")
    require("\\(" not in text and "\\[" not in text, "non-project math delimiter found")


def main() -> None:
    check_structure()
    text = manuscript_text()
    check_artifacts(text)
    check_citations(text)
    check_figures(text)
    check_style(text)
    print("PASS paper contract: structure, constants, conventional citations, figures, and prose gates")


if __name__ == "__main__":
    main()
