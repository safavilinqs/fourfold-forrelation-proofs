# Source snapshot

## Provenance

The snapshot was copied from repository state `6aff0c9` (`Close Round 4 theorem and feasibility decision`) on 2026-07-19 and audited on 2026-08-15. The audit modifies the Round-4 outward ledger, adaptive scope, artifacts, and verification entry points to remove unproved unbalanced routing placeholders. Round-2, Round-3, and other pre-audit Round-4 documents are retained as provenance; `../AUDIT_CORRECTION.md` supersedes their broad cross-number-sector claims and old numerical total. The relative directory structure required by imports is preserved.

| source tree | files | approximate size |
|---|---:|---:|
| `open_problem_forr4_passive_floor_consolidation_round_2` | 85 | 0.4 MB |
| `open_problem_forr4_passive_floor_consolidation_round_3` | 136 | 1.4 MB |
| `open_problem_forr4_passive_floor_consolidation_round_4` | 248 | 2.8 MB |
| complete version-controlled snapshot | 469 | 4.5 MB |

Python bytecode directories and local caches are excluded from these counts and from version control because they are generated files rather than theorem inputs.

## Principal artifact hashes

| artifact | SHA-256 |
|---|---|
| complete outward ledger | `33555d6aa5360971dd2783eec15b6929cf6d31e3f037315d6d647828b5213d44` |
| adaptive tree-frontier certificate | `41e19f7056673a1883c82d868ef085669cd82068a62e82947c571d678c5a13ca` |
| active six-dose resource certificate | `d592085c097cfd66ec761cb182cd0de3a9218ac54b4dc7725596993d471bcf31` |

The filenames and claim-to-artifact mapping are given in `../dev/SOURCE_MAP.md` and Appendix E of the manuscript. These hashes record the corrected 2026-08-15 focused run. The verification commands consume the included artifacts and source files directly; no network access or external data download is required.
