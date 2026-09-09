#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
PROJECT="$(cd "$ROOT/.." && pwd)"
ROUND4="$ROOT/source_snapshot/open_problem_forr4_passive_floor_consolidation_round_4"
PYTHON_BIN="${PYTHON_BIN:-python3}"
if [[ $# -gt 1 || ( $# -eq 1 && "$1" != "--full" ) ]]; then
  echo "usage: ./run_all.sh [--full]" >&2
  exit 2
fi
export PYTHONDONTWRITEBYTECODE=1
"$PYTHON_BIN" "$ROOT/check_paper_contract.py"
"$PYTHON_BIN" "$ROUND4/tests/active_six_resource_row.py"
"$PYTHON_BIN" "$ROUND4/tests/q64_complete_outward_ledger.py"
if [[ $# -eq 1 ]]; then
  echo "ARCHIVAL replay: old adaptive verdicts are superseded by AUDIT.md."
  PYTHON_BIN="$PYTHON_BIN" "$ROUND4/run_round4_checks.sh"
fi
echo "PASS current parallel-probe verification"
