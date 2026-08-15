VERDICT: CERTIFIED

# Independent audit of the twelve dual-endpoint entries

Date: 2026-07-17

The twelve entries in `q64_dual_endpoint_schur_insertion.py` have a valid arbitrary-correlated-diagonal one-batch coefficient. The original note omitted the essential factorization of the completed link touching a split support; this audit supplies that factorization and independently checks the physical endpoint slices, masks, reversal map, and link moments.

## Physical matrices

Write each physical support as $U_j=A_j\mathbin{\dot\cup}B_j$, where $A_j$ is in the occurrence-matrix row, $B_j$ is in its column, and $A_j\cap B_j=\varnothing$. A link moment is

$$
M_{a,b}(U,V)=\mathbb E_\omega L_U(\omega)R_V(\omega),
$$

where $L_U=\prod_{u\in U}L_u$ and $R_V=\prod_{v\in V}R_v$ are unit-modulus signed-permutation monomials. The occurrence entry is the product of the three consecutive link moments. Row and column probability laws may be arbitrary and correlated within their complete configurations.

The twelve entries form the complement/reversal closures of the following three representatives.

### Separate favorable singletons

For $(3,1,1,5):(1,0,1,3)$, rows are $(x,c,F)$ and columns are $(E,b,G)$ with $|E|=|G|=2$, $|F|=3$, $x\notin E$, and $F\cap G=\varnothing$. The physical matrix is

$$
K_{(x,c,F),(E,b,G)}=M_{3,1}(\{x\}\mathbin{\dot\cup}E,b)M_{1,1}(b,c)M_{1,5}(c,F\mathbin{\dot\cup}G).
$$

The cubic mask $x\notin E$ belongs to the first endpoint factor, the quintic mask $F\cap G=\varnothing$ belongs to the last endpoint factor, and the completed middle link is

$$
M_{1,1}(b,c)=\langle L_b,R_c\rangle.
$$

### Shared column singleton

For $(3,1,5,3):(1,0,2,3)$, rows are $(x,F,D)$ and columns are $(E,b,G)$ with $|E|=|F|=2$, $|G|=|D|=3$, $x\notin E$, and $F\cap G=\varnothing$. The physical matrix is

$$
K_{(x,F,D),(E,b,G)}=M_{3,1}(\{x\}\mathbin{\dot\cup}E,b)M_{1,5}(b,F\mathbin{\dot\cup}G)M_{5,3}(F\mathbin{\dot\cup}G,D).
$$

The last link is not a cross Gram merely because it is a moment. Its valid completed-link factorization is instead

$$
M_{5,3}(F\cup G,D)=\mathbb E[L_FL_GR_D]=\langle L_FR_D,L_G\rangle.
$$

The row feature depends on the complete row $(x,F,D)$ and the column feature depends on the complete column $(E,b,G)$. Both have norm one. No distinctness mask is assigned to this completed factor; the full quintic mask is already present in $M_{1,5}(b,F\mathbin{\dot\cup}G)$.

### Shared row singleton

For $(3,3,1,5):(0,2,1,3)$, rows are $(A,c,F)$ and columns are $(D,x,G)$ with $|A|=|G|=2$, $|D|=|F|=3$, $x\notin A$, and $F\cap G=\varnothing$. The physical matrix is

$$
K_{(A,c,F),(D,x,G)}=M_{3,3}(D,A\mathbin{\dot\cup}\{x\})M_{3,1}(A\mathbin{\dot\cup}\{x\},c)M_{1,5}(c,F\mathbin{\dot\cup}G).
$$

Here the completed first link has the exact factorization

$$
M_{3,3}(D,A\cup\{x\})=\mathbb E[L_DR_AR_x]=\langle R_A,L_DR_x\rangle.
$$

Again the cubic mask is carried by $M_{3,1}$ and the quintic mask by $M_{1,5}$, while the completed factor has unit feature norms.

## Endpoint factors

For every fixed cubic pair $E$ and adjacent singleton $b$, direct signed-permutation counting gives

$$
\sum_{x\notin E}|M_{3,1}(\{x\}\mathbin{\dot\cup}E,b)|^2\le E_2,
\qquad
E_2=\frac{q^2-2q+2}{q^2(q-1)}.
$$

Using the standard basis on the varying $x$ side and the displayed physical column vector on the fixed $(E,b)$ side gives a Schur feature factor at most $\sqrt{E_2}$. This is a complete physical slice and includes $x\notin E$.

For every fixed quintic triple $F$ and adjacent singleton $c$, the physical complementary-pair slice is

$$
\sum_{G\cap F=\varnothing,\ |G|=2}|M_{1,5}(c,F\mathbin{\dot\cup}G)|^2\le F_3,
$$

with

$$
F_3=\binom{q-3}{2}\frac1{q^2}+q(q-1)\left(\frac1{q^2}+\frac{q-4}{q^2(q-1)^2}\right)+(q-1)\binom q2\frac1{q^2(q-1)^2}.
$$

The same standard-basis factorization gives a Schur feature factor at most $\sqrt{F_3}$ and includes the full mask $F\cap G=\varnothing$.

If the two endpoint factors share the singleton $b$ or $c$, their Schur product uses tensor-product features at the same complete row or column index. For factorizations $A_{rc}=\langle u_r,v_c\rangle$ and $B_{rc}=\langle x_r,y_c\rangle$,

$$
(A\circ B)_{rc}=\langle u_r\otimes x_r,v_c\otimes y_c\rangle.
$$

Reusing a coordinate therefore does not duplicate a marginal law or assume independence. The arbitrary laws are applied once, after all three factors have been composed.

Combining the two endpoint factors and the appropriate unit completed-link factor gives

$$
\gamma_2(K)\le\sqrt{E_2F_3}.
$$

For arbitrary row and column probability laws $p,w$, any such feature factorization gives

$$
\left\|\operatorname{diag}(\sqrt p)K\operatorname{diag}(\sqrt w)\right\|_1\le\left(\max_r\|u_r\|\right)\left(\max_c\|v_c\|\right),
$$

by the Frobenius factorization bound. No product, uniform, exchangeable, or translation-invariant law is used.

Complementing swaps the occurrence row and column and hence transposes the weighted matrix. Reversal also transposes every physical cell coordinate, under which $M_{a,b}(U,V)=M_{b,a}(V^\top,U^\top)$. Both operations preserve the two simple-support masks and the trace norm, so the three representatives certify all twelve entries.

## Exact $q=64$ coefficient

At $q=64$,

$$
E_2=\frac{1985}{129024},
\qquad
F_3=\frac{62527}{43008},
$$

and therefore

$$
E_2F_3=\frac{124116095}{5549064192}.
$$

The committed coefficient is the first binary64 number whose exact square dominates this rational:

$$
c_{\rm dual}=0.14955611574342903<1.
$$

## Independent validation

`tests/q64_dual_endpoint_independent_audit.py` does not obtain its expected endpoint energies from the production formulas. It constructs all $384$ signed permutations at $q=4$, directly enumerates every simple-support fixed-pair cubic slice and fixed-triple quintic slice, and obtains $E_2=5/24$ and $F_3=7/8$. Its largest stored physical link matrix has $2{,}446{,}080$ signed-integer entries, below $5$ MB, and it constructs no giant occurrence matrix.

The test also checks 192 exact completed-link Gram identities, 120 physical reversal identities with coordinate transposition, direct $q=8$ permutation sums for all five link types, all twelve complement/reversal matrices, and adversarial nonuniform diagonal laws on sparse physical submatrices. These numerical checks do not replace the feature proof, but they independently protect every identity and mask placement used by it.

## Scope

This verdict removes the twelve-entry caveat and yields 888 of 888 supported balanced entries for the one-batch coefficient registry. It does not prove adaptivity. The next completed task is the dependency-exact outward-rounded $q=64$ ledger; adaptivity is eligible only if that certified total is below $1/3-10^{-3}$.
