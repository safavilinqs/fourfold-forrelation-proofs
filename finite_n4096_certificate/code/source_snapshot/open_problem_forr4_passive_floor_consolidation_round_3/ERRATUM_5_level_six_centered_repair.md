# Erratum 5 of record — level-six centered repair (2026-07-17)

*Filed by the coordinator from the independent rederivation program (clean-room Team Alpha, note R-106 with falsifier F-64/F-66 context; stop-the-line honored). This is the fifth genuine defect found in the internal corpus's history and the first found by the external rederivation rather than internal consolidation.*

**Defect.** In the level-six tree centered repair (`notes/low_level_terminal_centered_repair.md` and `searches/low_level_terminal_centered_repair.py`), the leftward centered expansion with fresh-$a_1$ cancellation has **no valid Stein instance**: the forced $\gamma'$ lies on the $b_0$ boundary-2 side (witnessed at the code's own site `('A',(1,0))`, by fresh enumeration, and by quadrature), independently of the layer-1 latents. Additionally the dangerous set has 33 elements, not 31 (a rank filter incorrectly drops two rank-$\le1$ partitions). The claim "exact re-expansion gives $N^{-1}$ on every retained branch" is therefore FALSE for the level-six tree as written.

**Repair (proved in the rederivation corpus, team_alpha R-106/R-108).** The corrected rightward expansion is valid and leaves 26 level-eight cases at decay $N^{-1/2}$ (rather than $N^{-1}$).

**Consequence for the theorem: none to the exponent.** With the repair, level six binds at $(1/2)/6=1/12$, jointly with the already-binding level twelve; the passive floor remains $D_{\mathsf P}=\Omega(N^{1/12})$ (rederivation note R-108, "S26′"). The 1/24 fallback is not engaged. Intermediate statements quoting per-branch $N^{-1}$ for the level-six tree must cite this erratum.

**Status of the other centered repairs.** The level-ten forest, both level-nine trees, and the level-seven tree were verified as written by the same referee pass (200/200 histories; dangerous sets complete; minima confirmed). The terminal-image enumeration (236/39/22, initials (12,8,8,4)) was independently re-implemented and reproduced exactly.
