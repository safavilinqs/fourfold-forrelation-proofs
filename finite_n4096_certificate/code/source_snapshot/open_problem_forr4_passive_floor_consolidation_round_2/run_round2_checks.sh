#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

here="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
python_bin="${PYTHON_BIN:-python3}"

if ! "$python_bin" -c 'import numpy, scipy, sympy' >/dev/null 2>&1; then
  recorded_python="/opt/homebrew/Caskroom/miniconda/base/bin/python3"
  if [[ -x "$recorded_python" ]]; then
    python_bin="$recorded_python"
  else
    printf 'NumPy, SciPy, or SymPy is unavailable under %s; set PYTHON_BIN to the audit interpreter.\n' "$python_bin" >&2
    exit 2
  fi
fi

printf 'python %s\n' "$($python_bin --version 2>&1)"
"$python_bin" "$here/tests/fiberwise_bilateral_stress.py"
"$python_bin" "$here/tests/collision_aware_stress.py"
"$python_bin" "$here/tests/insertion_collision_ledger.py"
"$python_bin" "$here/tests/reverse_tree_operator_stress.py"
"$python_bin" "$here/tests/open_frontier_operator_stress.py"
"$python_bin" "$here/tests/adaptive_mark_assignment_ledger.py"
"$python_bin" "$here/tests/temporal_square_function_barrier.py"
"$python_bin" "$here/tests/joint_probe_square_mass.py"
"$python_bin" "$here/tests/joint_occurrence_profile_packing.py"
"$python_bin" "$here/tests/two_copy_square_function_counterexample.py"
"$python_bin" "$here/tests/two_batch_trace_norm.py"
"$python_bin" "$here/tests/quadratic_bent_collision_barrier.py"
"$python_bin" "$here/tests/quadratic_endpoint_weighted_bound.py"
"$python_bin" "$here/tests/weighted_link_gram_contraction.py"
"$python_bin" "$here/tests/three_link_weighted_path_contraction.py"
"$python_bin" "$here/tests/gram_dressed_tail_contraction.py"
"$python_bin" "$here/tests/attenuated_exact_plant_variance.py"
"$python_bin" "$here/tests/two_node_chain_physical_stress.py"
"$python_bin" "$here/tests/mixed_component_physical_stress.py"
"$python_bin" "$here/tests/singleton_ancilla_norm_gap.py"
"$python_bin" "$here/tests/global_assignment_dichotomy_check.py"
"$python_bin" "$here/tests/all_singleton_physical_stress.py"
"$python_bin" "$here/tests/all_singleton_masked_graph_norm.py"
"$python_bin" "$here/tests/mixed_component_n4_spotcheck.py"
"$python_bin" "$here/searches/permutation_match_budget.py"
"$python_bin" "$here/searches/minimal_chain_recording_counterexample.py"
"$python_bin" "$here/searches/permutation_sector_spectra.py"
"$python_bin" "$here/searches/signed_permutation_full_sector_spectra.py"
"$python_bin" "$here/searches/signed_permutation_cubic_schur_lift.py"
"$python_bin" "$here/searches/signed_permutation_middle_cubic_schur_lift.py"
"$python_bin" "$here/searches/dose_six_sector_frontier.py"
"$python_bin" "$here/searches/signed_permutation_quintic_middle_bound.py"
"$python_bin" "$here/searches/signed_permutation_double_cubic_entries.py"
"$python_bin" "$here/searches/double_endpoint_slice_energies.py"
"$python_bin" "$here/searches/adjacent_double_cubic_slice_energies.py"
"$python_bin" "$here/searches/adjacent_double_cubic_orbit_scaling.py"
"$python_bin" "$here/searches/fixed_split_occupation_barrier.py"
"$python_bin" "$here/searches/double_endpoint_joint_schur_benchmark.py"
"$python_bin" "$here/tests/one_batch_minimal_chain_certificate.py"
printf 'PASS round-two reverse-tree stress checks\n'
