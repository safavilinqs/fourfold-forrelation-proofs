# Audit correction: total-signal-number scope

Date: 2026-08-15.

## Finding

The source state copied from commit `6aff0c9` claimed a passive lower bound for fresh batches with arbitrary photon-number coherence. That scope was not supported by its ledger.

For profile--split compatibility

$$
m=n+a-2t,
$$

the condition $2|t|=|a|$ follows only when the ket and bra occupations have equal total signal photon number. The original at-most-six state space also contains different-total-number pairs. An exact audit found 272 unbalanced high-sector profile-split/state incidences, forming 136 undirected occupation edges. The earlier code evaluated those edges using two-tier routing placeholders whose own source explicitly stated that they were not theorems.

Consequently, a passing original regression did not certify the advertised vacuum--nonvacuum-coherent passive class.

## Corrected certified statement

The finite theorem now applies to classically adaptive passive trees for which every fresh signal--idler state obeys

$$
[\rho,\widehat N_{\rm sig}\otimes I_{\rm id}]=0.
$$

The class permits arbitrary mixtures of total signal photon number, arbitrary idlers, repeated modes, arbitrary entanglement and coherence within a fixed-number sector, collective measurements, rare outcomes, branch-dependent dose partitions, and unrestricted classical feed-forward. It excludes coherence between different total signal-number sectors and coherent quantum memory between batches.

For this class, the theorem coefficient map contains exactly the 888 balanced high-sector entries and rejects every unbalanced routing value. Directed reconstruction gives

$$
P_{\rm Perron}\le
0.2587440963855772226792307915292556241,
$$

and, after the unchanged promise-conditioning cost,

$$
\operatorname{TV}(P_+,P_-)
\le
0.2609692247922079249341809573938165614.
$$

The adaptive direct-sum factorization has multiplier one within the corrected class.

## Open extension

The same lower bound for fresh passive probes with coherence between distinct total signal-number sectors is open. Restoring that claim requires theorem-backed physical coefficients for every unbalanced occurrence kernel that can connect the at-most-six occupation states. Passing inventory, artifact-consistency, and random stress tests cannot substitute for those bounds.

Historical Round-2, Round-3, and pre-audit Round-4 planning or run documents under `code/source_snapshot/` are provenance. If they state the broader theorem or the old numerical total, this correction supersedes them.
