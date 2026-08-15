# Paper-closeout theorem package

STATUS: NUMBER-SECTOR-INCOHERENT MATHEMATICAL SEPARATION CERTIFIED AT $N=4096$

EXPERIMENTAL LABEL: IDEAL RESOURCE SPECIFICATION; PRESENT-DAY PLATFORM CREDIBILITY NOT CERTIFIED

Date: 2026-07-18

## Task

Let $N=4096=64^2$, let $H=H_N$ be the normalized Sylvester matrix, and let $x^{(1)},\ldots,x^{(4)}\in\{\pm1\}^N$. Define

$$
F_{4,H}(x)=\frac1N(x^{(1)})^THD_2HD_3Hx^{(4)},
\qquad
D_j=\operatorname{diag}(x^{(j)}).
$$

The promise problem distinguishes

$$
F_{4,H}(x)\ge\frac14
\qquad\text{from}\qquad
F_{4,H}(x)\le-\frac14.
$$

There are four sign banks and therefore $M=4N=16{,}384$ unknown sign coordinates.

## Hard pair and conditioning

Write $q=64$ and let $S_q=\sqrt q\,H_q$ be the unnormalized Sylvester sign matrix. A uniform signed permutation is specified by a uniform permutation $\pi$ of $[q]$ and independent uniform signs $\sigma_b$. Associate to it two $q\times q$ sign arrays

$$
L_{a,b}=S_q(a,\pi(b))\sigma_b,
\qquad
R_{\pi(b),a}=\sigma_bS_q(b,a),
$$

and flatten each array into a length-$N$ sign string. Draw three independent signed permutations $P_1,P_2,P_3$ and form the exact positive plant

$$
z^{(1)}=L(P_1),
\quad
z^{(2)}=R(P_1)\odot L(P_2),
\quad
z^{(3)}=R(P_2)\odot L(P_3),
\quad
z^{(4)}=R(P_3).
$$

This plant satisfies $F_{4,H}(z)=1$. The exact negative plant is obtained by replacing $z^{(1)}$ by $-z^{(1)}$ and satisfies $F_{4,H}(z)=-1$.

Fix $\beta=19/25$. Independently for every coordinate in every block, multiply the exact plant by a sign $\eta$ with

$$
\Pr(\eta=1)=\frac{1+\beta}{2},
\qquad
\Pr(\eta=-1)=\frac{1-\beta}{2}.
$$

The positive hard law is this attenuated plant conditioned on $F_{4,H}\ge1/4$; the negative hard law is its reflected counterpart conditioned on $F_{4,H}\le-1/4$. The outward upper on the sum of the two bad-promise probabilities is

$$
0.0022251284066307022549501658645609374.
$$

## Passive theorem

Consider every finite classically adaptive passive tree that, at each node, prepares a fresh signal--idler state block diagonal in total signal photon number; passes its signal through the diagonal sign oracle; makes an arbitrary collective batch POVM; and chooses the next fresh batch from the complete classical outcome history. Within each fixed-number sector, signal entanglement, idlers, repeated modes, and coherence are arbitrary. Charge every photon pass through an unknown sign coordinate and require total dose at most six on every root-to-leaf branch. Coherence between different total signal-number sectors and coherent quantum memory or control between batches are excluded.

For the two conditioned hard laws above, every such passive strategy obeys

$$
\operatorname{TV}(P_+,P_-)
\le
0.2609692247922079249341809573938165614
<
\frac13-10^{-3}.
$$

The certified reserve to the acceptance threshold is

$$
0.0713641085411254083991523759395167720.
$$

With equal priors, the average error is therefore at least

$$
\frac{1-\operatorname{TV}(P_+,P_-)}2
\ge
0.3695153876038960375329095213030917193
>
\frac13.
$$

Consequently some promised input has error greater than $1/3$ for every passive hard-dose-six strategy, so

$$
D_{\mathsf P}^{\rm hard}(4096)>6
$$

in the declared classically adaptive model.

## Active theorem

Let $|u\rangle=N^{-1/2}\sum_i|i\rangle$ and define

$$
|L_x\rangle=D_2HD_1|u\rangle,
\qquad
|R_x\rangle=HD_3HD_4|u\rangle.
$$

One photon coherently follows the chronological branch words $D_1,H,D_2$ and $D_4,H,D_3,H$. A path-$X$ measurement gives a binary flag with expectation $F_{4,H}(x)$. Three independent flags followed by majority vote use deterministic branchwise hard dose six, retain every outcome, and have worst promised error

$$
\frac{81}{256}=0.31640625<\frac13,
$$

with exact margin $13/768$. Hence

$$
D_{\mathsf A}^{\rm hard}(4096)\le6.
$$

## Finite-size separation corollary

Under the same four-bank sign oracle and the same branchwise photon-pass meter,

$$
\boxed{
D_{\mathsf P}^{\rm hard}(4096)>6
\quad\text{and}\quad
D_{\mathsf A}^{\rm hard}(4096)\le6.}
$$

The passive statement permits unrestricted classical feed-forward between fresh quantum batches that are block diagonal in total signal photon number. The active statement is pointwise for every promised input. Neither statement discounts loss through postselection. The extension to fresh probes with coherence between different total signal-number sectors is open.

## Resource row and experimental boundary

| quantity | certified value |
|---|---:|
| problem dimension | $N=4096$ |
| unknown sign coordinates | $M=16{,}384$ |
| active photons | 3 sequential or parallel |
| charged sign traversals | 2 per photon; 6 total |
| public transforms | one $H_{4096}$ on the left branch; two on the right |
| receiver | mode-insensitive path-$X$ port detection and majority vote |
| postselection | none |
| active error | $81/256$ |
| passive transcript upper | $0.260969224792208$ outward |
| passive access boundary | number-sector-incoherent fresh batches; classical memory only between batches |

The active implementation requires combined zero-bias retained contrast strictly above $0.904294855157$. `EXPERIMENTAL_FEASIBILITY_DECISION.md` finds no reviewed platform demonstration of the required $4096$-dimensional coherent Sylvester transform at that threshold. The paper may claim a rigorous finite-size theory separation with an explicit ideal optical resource row. It should not call the $N=4096$ realization experimentally demonstrated or presently credible without new device evidence.

## Proof dependency map

| claim | authoritative dependency |
|---|---|
| all physical one-batch coefficients | `Q64_MASKED_UNIVERSAL_AUDIT.md` and its dependency-exact registry |
| outward occurrence and promise total | `Q64_COMPLETE_OUTWARD_LEDGER.md` and `artifacts/q64_complete_outward_ledger.json` |
| arbitrary classical adaptivity | `Q64_ADAPTIVE_TREE_FRONTIER_THEOREM.md` and `artifacts/q64_adaptive_tree_frontier.json` |
| active dose and error | `ACTIVE_SIX_DOSE_RESOURCE_ROW.md` and `artifacts/active_six_resource_row.json` |
| active contrast threshold | `ACTIVE_SIX_ROBUSTNESS_GATE.md` and `artifacts/active_six_robustness_gate.json` |
| hardware decision | `EXPERIMENTAL_FEASIBILITY_DECISION.md` and `artifacts/q64_experimental_feasibility_gate.json` |

## Reproduction

From the Round-4 folder run:

    ./run_round4_checks.sh

The command regenerates and compares the committed mathematical artifacts and runs the complete inherited and Round-4 regression suites.

## Nonclaims

This package does not prove a lower bound for fresh probes with coherence between total signal-number sectors or against protocols carrying coherent quantum memory between batches, does not prove optimal active or passive dose, does not certify a laboratory implementation, and does not make the signed-permutation hard distribution part of the task definition. The hard distribution is only the lower-bound witness. The audit found 272 unbalanced high-sector split/state incidences on 136 different-number occupation edges; their coefficients remain open.
