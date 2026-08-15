# Opposite cubic--quintic endpoint orbit factorization

Status: exact endpoint formulas and fixed-orbit contraction proved.  The
subsequent mixed-orbit witness in
`opposite_endpoint_mixed_orbit_obstruction.md` shows that this contraction
does not extend to arbitrary diagonal laws at the required coefficient.

## 1. Target

The first unresolved reversal pair is

$$
(3,1,1,5)/(5,1,1,3).
$$

For the critical split $(2,0,1,2)$, write the rows as $(A,C,c)$ and the
columns as $(e,D,b)$, where $A,C$ are selected endpoint pairs, $e$ is the
remaining cubic endpoint cell, $D$ is the remaining quintic endpoint triple,
and $b,c$ are the middle singletons.  The direct squared-slice estimate is
$0.780328871525$ and is unusable.

This note keeps the endpoint signs and the middle Walsh transform together.
It gives a much smaller exact answer after the pair differences and triple
translation shape are fixed.

## 2. Exact xor-labelled endpoint moments

Let $H_N$ be the normalized Walsh matrix and let $Q$ be an endpoint support.
The cubic and quintic endpoint moments have the form

$$
M_{31}(Q,z)=v_3(Q)H_N(\mathop{\rm xor}Q,z),
\qquad
M_{51}(Q,z)=v_5(Q)H_N(\mathop{\rm xor}Q,z).
\tag{2.1}
$$

For a cubic support, $v_3=1$ when all three cells use one hidden column,
$v_3=-1/(q-1)$ when they use two columns, and $v_3=0$ when they use three.

For a quintic support with one odd hidden column, take the row XOR in each
even column.  There can be zero, one, or two such values, denoted by
$\alpha,\beta$.  Then

$$
v_5=
\begin{cases}
1,&\text{no even column},\\
1,&\text{one even column and }\alpha=0,\\
-1/(q-1),&\text{one even column and }\alpha\ne0,\\
-1/(q-1),&\text{two even columns and }\alpha=\beta,\\
2/[(q-1)(q-2)],&\text{two even columns and }\alpha\ne\beta.
\end{cases}
\tag{2.2}
$$

The value is zero if the support has more than one odd hidden column.  These
formulas follow from the one- and two-character averages of an injective
signed-permutation column assignment.  Exhaustive $q=4$ enumeration checks
(2.1)--(2.2) entry by entry for every degree-three and degree-five support.

## 3. Fixed translation-orbit factorization

Fix the nonzero XORs $x$ and $y$ of the selected pairs $A$ and $C$, and fix
one translation orbit $\mathcal T$ for the triple $D$.  Define the endpoint
incidence matrices

$$
B_3[A,e]=v_3(A\cup\{e\}),
\qquad
B_5[C,D]=v_5(C\cup D),\qquad D\in\mathcal T,
\tag{3.1}
$$

with intersecting supports assigned value zero.  The physical block is

$$
K[(A,C,c),(e,D,b)]
=M_{31}(A\cup\{e\},b)H_N(b,c)M_{51}(C\cup D,c).
\tag{3.2}
$$

Because $\mathop{\rm xor}A=x$, $\mathop{\rm xor}C=y$, and the triple shape
is fixed, Walsh orthogonality in $b$ gives the exact row Gram

$$
KK^*= {I_N\over N^2}\otimes(B_3B_3^*)\otimes(B_5B_5^*).
\tag{3.3}
$$

The block has $N^3/4$ rows and $N^3$ columns.  Therefore its normalized
nuclear coefficient is

$$
\boxed{
{2\over N^3}
\operatorname{tr}\sqrt{B_3B_3^*}
\operatorname{tr}\sqrt{B_5B_5^*}.}
\tag{3.4}
$$

Each endpoint matrix is a pair-orbit-by-translation convolution.  If $f$ is
its response from the base pair, its singular values are
$|\widehat f(\xi)|/\sqrt2$ over the $N/2$ frequencies satisfying
$\langle\xi,x\rangle=0$.  Thus (3.4) is evaluated with length-$N$ Walsh
transforms; no $N^3$ matrix is formed.  A direct $q=4$ construction of the
$1024$-by-$4096$ physical matrix agrees with (3.3) through numerical
precision.

## 4. Aligned orbit value

For aligned vertical pairs and a triple consisting of a matching vertical
pair plus one odd cell, put $A=q(q-2)$.  The high-amplitude part alone has

$$
C_{\rm high}(q)={4A^2\over q^6},
\tag{4.1}
$$

and retaining all signed endpoint amplitudes gives

$$
\boxed{
C_{\rm full}(q)=
{4[A+2][A+2(q-2)/(q-1)]\over q^6}.}
\tag{4.2}
$$

The values are:

| $q$ | $N$ | $C_{\rm high}$ | $C_{\rm full}$ |
|---:|---:|---:|---:|
| 4 | 16 | 0.0625000000 | 0.0911458333 |
| 8 | 64 | 0.0351562500 | 0.0379289900 |
| 16 | 256 | 0.0119628906 | 0.0121702830 |
| 32 | 1024 | 0.00343322754 | 0.00344731635 |

Exhaustive scans over all pair differences and all triple translation shapes
at $q=4$ and $q=8$ find (4.1)--(4.2) as the largest *single fixed-orbit*
values.  This is evidence about the right extremizer, not yet a proof for
arbitrary $q$.

## 5. Finite-size significance and remaining gap

If both reversal profiles could be charged the $q=32$ value in (4.2), while
leaving every other open profile omitted, the compatible Perron diagnostic
would give

$$
F+2\epsilon_\beta=0.332197918529<1/3,
$$

with slack $0.001135414804$.  The largest uniform coefficient these two
profiles can tolerate is $0.006241223680$, so the fixed-orbit value has a
factor $1.81$ of headroom.

This is not an accepted ledger entry.  Translation twirling permits a law to
mix pair differences and triple shapes, and concavity means that the best
mixture need not be a pure orbit.  A triangle inequality over individual
orbits destroys the contraction.  The follow-up exact mixed-orbit reduction
finds an explicit $q=32$ coefficient $0.0395939553$, so the hoped-for
uniform upper below $0.00624122$ is false.

Before that witness was found, the possible next outcomes were:

1. give a supporting-dual upper below $0.00624122$ for the critical split;
2. give a small type matrix whose Perron or Collatz bound is below that
   threshold; or
3. produce a mixed-orbit witness above threshold, quantitatively closing
   this route and triggering the planned pivot — this is what occurred.

The independent scalar ledger now needs a promise or accepted-sector repair,
a joint cut/profile contraction, or a hard-instance pivot before this
reversal pair can enter a finite-size theorem.

## 6. Reproduction

- `searches/opposite_endpoint_orbit_scan.py` implements the exact endpoint
  amplitudes, length-$N$ response spectra, physical $q=4$ comparison, and
  $q=4,8,16,32$ aligned values.  Pass `--exhaustive` to reproduce the full
  $q=4,8$ pure-orbit scans.
- `tests/occupation_compatible_sector_optimization.py` exhaustively checks
  the endpoint formulas at $q=4$, checks (4.1)--(4.2), and scans every $q=4$
  fixed orbit.
- `./run_round3_checks.sh` runs the inherited and Round 3 regressions.
