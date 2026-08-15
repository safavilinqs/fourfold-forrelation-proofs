#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

here="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
round_two="$here/../open_problem_forr4_passive_floor_consolidation_round_2"

if [[ ! -x "$round_two/run_round2_checks.sh" ]]; then
  printf 'Missing executable round-two baseline: %s\n' "$round_two/run_round2_checks.sh" >&2
  exit 2
fi

"$round_two/run_round2_checks.sh"

shopt -s nullglob
round_three_tests=("$here"/tests/*.py)
if ((${#round_three_tests[@]})); then
  python_bin="${PYTHON_BIN:-/opt/homebrew/Caskroom/miniconda/base/bin/python3}"
  for test_file in "${round_three_tests[@]}"; do
    "$python_bin" "$test_file"
  done
fi

printf 'PASS round-three initialization and inherited baseline\n'
