#!/usr/bin/env python3
"""Check current manuscript bounds against the frozen finite certificate."""
import json
import re
from fractions import Fraction
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT.parent
ARTIFACTS = PROJECT / "code/source_snapshot/open_problem_forr4_passive_floor_consolidation_round_4/artifacts"


def main():
    text = (PROJECT / "main.tex").read_text()
    payload = json.loads((ARTIFACTS / "q64_complete_outward_ledger.json").read_text())
    ledger = payload["result"]
    active = json.loads((ARTIFACTS / "active_six_resource_row.json").read_text())["result"]
    displayed = {
        "total_upper": "0.260969224792207925",
        "collatz_perron_upper": "0.258744096385577223",
        "promise_loss_upper": "0.002225128406630703",
    }
    for key, value in displayed.items():
        assert value in text, (key, "display missing")
        assert Fraction(value) >= Fraction(ledger[key]), (key, "inward rounding")
    total = Fraction(ledger["total_upper"])
    assert Fraction(ledger["collatz_perron_upper"]) + Fraction(ledger["promise_loss_upper"]) <= total
    error = Fraction("0.369515387603896037")
    assert "0.369515387603896037" in text
    assert Fraction(1, 3) < error <= (1 - total) / 2
    assert total < Fraction(1, 3) - Fraction(1, 1000)
    assert active["majority_error_exact"] == "81/256"
    assert active["total_hard_dose"] == 6
    assert Fraction(81, 256) == Fraction(3, 8)**3 + 3 * Fraction(5, 8) * Fraction(3, 8)**2
    assert ledger["certified_balanced_high_sector_coefficients"] == 888
    assert ledger["supported_balanced_entries"] == 888
    assert ledger["open_balanced_entries"] == 0
    assert ledger["excluded_unbalanced_high_sector_entries"] == 5128
    assert ledger["excluded_unbalanced_high_sector_incidence_records"] == 272
    assert ledger["excluded_unbalanced_high_sector_undirected_edges"] == 136
    for row in payload["coefficient_registry"]:
        assert 2 * sum(row["split"]) == sum(row["profile"])
    assert len(payload["collatz_candidate"]) == 210
    assert all(Fraction(v) > 0 for v in payload["collatz_candidate"])

    # Old multiplier-one artifacts are not inputs to the supported claim.
    for source in (text, (ROOT / "asymptotic_single_pass_floor/main.tex").read_text()):
        assert "parallel" in source and "gap" in source
    assert "block diagonal in total signal photon number" in text
    assert (ROOT / "AUDIT.md").is_file()
    assert (PROJECT / "code/COEFFICIENTS.md").is_file()
    bib = (PROJECT / "references.bib").read_text()
    keys = set(re.findall(r"@\w+\{([^,]+),", bib))
    for group in re.findall(r"\\cite(?:p|t)?\{([^}]+)\}", text):
        assert set(group.split(",")) <= keys
    assert "\\input{sections/" not in text and "\\input{appendix/" not in text
    print("PASS proof contract: outward arithmetic, registry, citations and scope")


if __name__ == "__main__":
    main()
