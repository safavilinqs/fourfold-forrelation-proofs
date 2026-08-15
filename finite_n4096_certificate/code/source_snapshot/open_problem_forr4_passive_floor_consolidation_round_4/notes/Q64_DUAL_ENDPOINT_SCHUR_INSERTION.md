# Dual endpoint-slice Schur insertion at q=64

Date: 2026-07-16; independently audited 2026-07-17

Status: certified arbitrary-correlated-diagonal one-batch coefficients for twelve entries in three complement/reversal orbits. `Q64_DUAL_ENDPOINT_INDEPENDENT_AUDIT.md` reconstructs all three physical matrices, supplies the completed-link Gram factorization omitted here, and validates the result by exact $q=4$ and direct $q=8$ tests. In the current registry these are the final twelve formerly caveated entries, so the one-batch balanced count is 888 of 888. The adaptive lift remains separate.

## Result

Twelve remaining split-cubic/split-quintic entries have both of the following
properties:

1. the split cubic is adjacent to a singleton lying on the cubic-pair side;
2. the split quintic is adjacent to a singleton lying on the quintic-majority
   side.

They form one degree-ten orbit and two degree-twelve orbits. All have a
balanced $3|2$ quintic split. Their common arbitrary-law coefficient is

$$
\boxed{0.149556115743}.
$$

After insertion, the routing values are

$$
\begin{aligned}
\beta&=0.746087134907,\\
P_{\rm Perron}&=0.312576926665,\\
P_{\rm promise}&=0.0173066793099,\\
P_{\rm total}&=0.329883605975.
\end{aligned}
$$

The routing margin is $0.00344972735879$. Because the previous target on
these entries was smaller than the proved coefficient, honest insertion
spends $0.000635351840339$ of provisional margin. The remaining quintic
inventory is 152 entries: 104 extreme and 48 balanced.

## Schur-multiplier composition

For a matrix symbol $A$, let $\gamma_2(A)$ be the least product of the
maximum row- and column-feature norms over factorizations
$A_{rc}=\langle u_r,v_c\rangle$. Schur multiplication obeys

$$
\|A\circ X\|_1\le\gamma_2(A)\|X\|_1,
$$

and the norm is submultiplicative under entrywise products.

For an internally split support block, the physical distinctness mask must
be assigned to one adjacent link. A bare cross-Gram completion is not enough.
The favorable singleton placements do exactly this for both split blocks.

### Cubic endpoint

Fix the cubic pair and its adjacent singleton. The vector over the remaining
cubic cell has squared norm bounded by the exact fixed-pair slice

$$
E_2=0.0153847346230.
$$

This factorization includes the cubic within-block distinctness mask and
gives

$$
\gamma_2(A_{\rm cubic})\le\sqrt{E_2}
=0.124035215254.
$$

### Quintic endpoint

Fix the quintic majority triple and its adjacent singleton. The vector over
the complementary pair includes the quintic distinctness mask and has
squared norm

$$
F_3=1.45384579613.
$$

Therefore

$$
\gamma_2(A_{\rm quintic})\le\sqrt{F_3}
=1.20575528037.
$$

### Remaining link

After assigning the cubic and quintic masks to those endpoint factors, every other link can be completed across overlaps. For the $(3,1,5,3)$ representative, the nontrivial identity is

$$
M_{5,3}(F\cup G,D)=\mathbb E[L_FL_GR_D]=\langle L_FR_D,L_G\rangle.
$$

For the $(3,3,1,5)$ representative it is

$$
M_{3,3}(D,A\cup\{x\})=\mathbb E[L_DR_AR_x]=\langle R_A,L_DR_x\rangle.
$$

These are genuine row/column feature factorizations on the complete occurrence indices. They remain valid across overlaps because the physical distinctness masks have already been retained in the endpoint factors. The third representative has the ordinary unit factor $M_{1,1}(b,c)=\langle L_b,R_c\rangle$.

Submultiplicativity now gives

$$
\gamma_2(K)
\le\sqrt{E_2F_3}
=0.149556115743.
$$

The argument also allows the two favorable endpoint factors to share the same singleton block. Pointwise Schur composition tensors the two feature vectors at that same complete row or column index; no marginal law is used twice and no independence assumption is introduced.

## Scope and next topology

Historically this theorem converted twelve routing targets into arbitrary-law bounds without changing the all-quintic two-tier proxy. Subsequent theorems have now resolved the other entries. The exact dual coefficient is derived from $E_2F_3=124116095/5549064192$ and rounded outward. The independent audit is the authoritative proof supplement.

Reproduce with:

    /opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/q64_dual_endpoint_schur_insertion.py --write-artifact
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_dual_endpoint_schur_insertion.py
    /opt/homebrew/Caskroom/miniconda/base/bin/python3 tests/q64_dual_endpoint_independent_audit.py
