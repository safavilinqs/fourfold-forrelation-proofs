#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
PROJECT="$(cd "$ROOT/.." && pwd)"
ROUND4="$ROOT/source_snapshot/open_problem_forr4_passive_floor_consolidation_round_4"

if [[ $# -gt 1 || ( $# -eq 1 && "$1" != "--full" ) ]]; then
  echo "usage: ./run_all.sh [--full]" >&2
  exit 2
fi

CACHE_ROOT="${TMPDIR:-/tmp}/finite_n4096_verification_cache"
mkdir -p "$CACHE_ROOT/matplotlib" "$CACHE_ROOT/xdg"
export MPLCONFIGDIR="$CACHE_ROOT/matplotlib"
export XDG_CACHE_HOME="$CACHE_ROOT/xdg"

if [[ -n "${PYTHON_BIN:-}" ]]; then
  PYTHON_BIN="$PYTHON_BIN"
else
  for candidate in \
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 \
    /opt/homebrew/opt/python@3.13/bin/python3.13 \
    python3; do
    if "$candidate" -c 'import numpy, scipy, sympy' >/dev/null 2>&1; then
      PYTHON_BIN="$candidate"
      break
    fi
  done
fi

if [[ -z "${PYTHON_BIN:-}" ]]; then
  echo "No Python interpreter with NumPy, SciPy, and SymPy was found; set PYTHON_BIN explicitly." >&2
  exit 2
fi

echo "Python interpreter: $PYTHON_BIN"
"$PYTHON_BIN" --version

PYTHONDONTWRITEBYTECODE=1 "$PYTHON_BIN" "$ROOT/check_paper_contract.py"
PYTHONDONTWRITEBYTECODE=1 "$PYTHON_BIN" "$PROJECT/figures/generate_figures.py" --check

FOCUSED_TESTS=(
  active_six_resource_row.py
  active_six_robustness_gate.py
  q64_complete_outward_ledger.py
  q64_adaptive_tree_frontier.py
  q64_experimental_feasibility_gate.py
  q64_paper_closeout_package.py
)

for test_name in "${FOCUSED_TESTS[@]}"; do
  PYTHONDONTWRITEBYTECODE=1 "$PYTHON_BIN" "$ROUND4/tests/$test_name"
done

if [[ "${1:-}" == "--full" ]]; then
  PYTHON_BIN="$PYTHON_BIN" "$ROUND4/run_round4_checks.sh"
fi

echo "PASS number-sector-incoherent paper verification: $PROJECT"
