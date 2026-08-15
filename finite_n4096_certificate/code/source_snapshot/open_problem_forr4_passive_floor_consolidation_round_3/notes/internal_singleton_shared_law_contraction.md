# Internal-singleton shared-law contraction

Date: 2026-07-15

Status: rigorous arbitrary-diagonal bound for the fourth balanced orbit.
The coefficient \(0.02509674612\) is below both the provisional \(1/32\)
and the route gate \(0.04504054678\).  The proof is the first Round 3
contraction that uses an exact shared-law endpoint coefficient rather than
multiplying two worst endpoint slices.

## 1. Target and the failed separate-slice bound

For profile \((3,1,1,5)\) and split \((1,1,1,2)\), write the cubic support
as \(\{x\}\dot\cup E\), \(|E|=2\), and the quintic support as
\(F\dot\cup G\), \(|F|=2\), \(|G|=3\).  Rows are \((x,b,c,F)\), columns are
\((E,G)\), and the exact tensor is

$$
K_{(x,b,c,F),(E,G)}
=M_{3,1}(\{x\}\cup E,b)H_N(b,c)M_{1,5}(c,F\cup G).
\tag{1.1}
$$

The singleton labels \(b,c\) are both on the row side, so
\(|H_N(b,c)|=1/q\) is an internal factor.

The direct symmetric row/column slice theorem gives

$$
c_{\rm slice}
=\sqrt{{E_1F_2\over N}}
=0.09914702581,
\tag{1.2}
$$

where

$$
E_1={N+2\over2N},
\qquad
F_2={159457\over7936}.
$$

This is rigorous but fails the finite-size gate.  Inserting (1.2) and
reoptimizing gives

$$
\operatorname{TV}_{\rm diagnostic}=0.33719291910,
\tag{1.3}
$$

which exceeds \(1/3\) by \(0.00385958577\).  A theorem that merely pairs
the two worst endpoint slices cannot preserve the route.

## 2. Exact cubic arbitrary-law coefficient

Define the split cubic endpoint matrix

$$
A_{(x,b),E}
=\mathbf 1_{\{x\notin E\}}
M_{3,1}(\{x\}\cup E,b).
\tag{2.1}
$$

For probability laws \(p,r\), write

$$
\Phi_A(p,r)
=\|D_p^{1/2}AD_r^{1/2}\|_1.
\tag{2.2}
$$

The function \(\Phi_A\) is jointly concave in \((p,r)\).  One way to see
this is the root-fidelity variational formula

$$
\operatorname{Tr}\sqrt{P^{1/2}QP^{1/2}}
=\frac12\inf_{X>0}
\left\{\operatorname{Tr}(PX)+\operatorname{Tr}(QX^{-1})\right\},
\tag{2.3}
$$

with \(P=D_p\) and \(Q=AD_rA^*\).  For fixed \(X\), the expression is
linear in both laws, so the infimum is jointly concave.

The affine hidden-row/hidden-column symmetries preserve \(A\) up to row and
column phases.  Twirling (2.2) therefore cannot decrease it.  The row law
becomes uniform, while the pair law retains only three weights:

- vertical pair difference;
- horizontal pair difference; and
- general pair difference.

For a representative difference \(d\), let \(\rho_d\) be the exact cubic
response and take its Walsh transform on the \(d\)-even frequencies.
At \(q=32\), after multiplying \(\rho_d\) by \(q-1\), the exact absolute
spectrum sums are

$$
S_{\rm v}=59644,\qquad
S_{\rm h}=1984,\qquad
S_{\rm g}=1984.
\tag{2.4}
$$

There are \(31,31,961\) differences of the three types.  Orthogonality in
\(b\) separates their Walsh blocks.  The three coefficients multiplying
the square roots of the orbit masses are

$$
\begin{aligned}
a_{\rm v}&={59644\sqrt{31}\over31\cdot32^3}
          =0.326915851870,\\
a_{\rm h}&={1984\sqrt{31}\over31\cdot32^3}
          =0.0108745397712,\\
a_{\rm g}&={1984\cdot31\over31\cdot32^3}
          =0.060546875.
\end{aligned}
\tag{2.5}
$$

Cauchy--Schwarz optimizes the three remaining orbit weights and proves the
exact arbitrary-law coefficient

$$
\boxed{
c_A=\sqrt{a_{\rm v}^2+a_{\rm h}^2+a_{\rm g}^2}
=0.332653203639.
}
\tag{2.6}
$$

Complete \(q=4\) construction attains the analogous twirled value
\(0.693158531651\), and random nonsymmetric diagonal laws remain below it.
This protects both the symmetry reduction and the direction of the
concavity argument.

## 3. Complete the quintic and subtract overlaps

Let

$$
B_{(F,c),G}
=\mathbf 1_{\{F\cap G=\varnothing\}}
M_{1,5}(c,F\cup G).
\tag{3.1}
$$

Complete it to every pair/triple pair by

$$
\widetilde B_{(F,c),G}
=\mathbb E[(X_cY_F)Y_G].
\tag{3.2}
$$

This is a unit-feature cross Gram kernel, so its trace-class Schur
multiplier norm is at most one.  On disjoint \(F,G\), it equals (3.1).
When they overlap, the repeated \(Y\)-coordinates cancel and (3.2) becomes
a cubic or singleton endpoint moment.  Hence

$$
B=\widetilde B-C,
\tag{3.3}
$$

where \(C\) is supported on \(F\cap G\ne\varnothing\).

Fix \(F=\{f_1,f_2\}\) and \(c\).  If \(|F\cap G|=1\), the two possible
shared elements leave a cubic support containing the other fixed element.
Their total squared energy is at most \(2E_1\).  If
\(|F\cap G|=2\), the remaining singleton contributes

$$
\sum_{g\notin F}|H_N(c,g)|^2=1-{2\over N}.
$$

Therefore every row of \(C\) has squared energy at most

$$
2E_1+1-{2\over N}
={N+2\over N}+1-{2\over N}
=2.
\tag{3.4}
$$

Using the row vector itself as a Schur feature proves

$$
\|C\|_{\rm Schur}\le\sqrt2,
\qquad
\boxed{\|B\|_{\rm Schur}\le1+\sqrt2.}
\tag{3.5}
$$

Complete \(q=4\) enumeration gives maximum collision energy
\(23/12<2\) and checks (3.3) entry by entry.

## 4. Shared-law whole-chain theorem

Repeat \(A\) over the ancillary labels \((c,F,G)\).  For arbitrary correlated
row and column laws, duplicate-row and duplicate-column compression
aggregates them to laws on \((x,b)\) and \(E\).  Equation (2.6) still bounds
the resulting weighted nuclear norm by \(c_A\).

Now Schur multiply by the repeated quintic kernel \(B\).  Equation (3.5)
costs at most \(1+\sqrt2\), without duplicating either physical law.  Finally
remove the row-only phase of \(H_N(b,c)\), retaining its modulus \(1/q\).
Thus

$$
\boxed{
\|D_p^{1/2}KD_r^{1/2}\|_1
\le {c_A(1+\sqrt2)\over q}
\sqrt{\left(\sum p\right)\left(\sum r\right)}.
}
\tag{4.1}
$$

At \(N=1024\),

$$
\boxed{
c_{(3,1,1,5):(1,1,1,2)}
=0.0250967461185
< {1\over32}
<0.0450405467778.
}
\tag{4.2}
$$

Complement and path reversal prove the same value on the other three cuts.

The mechanism matters: the cubic endpoint is optimized against one shared
law before the quintic is applied as a Schur dressing.  Multiplying the
separate worst slice factors would give the failing value (1.2).

## 5. Physical-law stress tests

The exact vertical-pair/vertical-triple mixture that is analogous to the
earlier opposite-endpoint witness gives

$$
c_{\rm vertical}=0.00665619197928
\tag{5.1}
$$

at \(q=32\).  Fixed translation orbits are smaller still; an aligned
vertical family gives \(0.000451712699\).  These are lower diagnostics, not
upper theorems, but neither approaches (4.2) or the ledger gate.  The
arbitrary-law proof, rather than the searches, establishes (4.2).

## 6. Ledger consequence

Insert (4.2) together with the three earlier balanced theorems.  Reoptimizing
the 210-state Perron diagnostic and the extended promise loss gives

$$
\beta=0.779338838167,\qquad
\operatorname{TV}_{\rm diagnostic}=0.332368300150.
\tag{6.1}
$$

The remaining slack is

$$
\boxed{
{1\over3}-\operatorname{TV}_{\rm diagnostic}
=0.000965033183.
}
\tag{6.2}
$$

Sixteen of 888 balanced entries now have chain-aware arbitrary-law
coefficients.  The largest remaining provisional orbit is

$$
(3,1,1,5):(0,0,1,4),
\tag{6.3}
$$

with Perron contribution \(0.00128192174\).  Reoptimizing while varying only
that four-cut orbit gives its next acceptance gate:

$$
\boxed{
c_{(3,1,1,5):(0,0,1,4)}
<0.0542506297760
=1.73602015283/q.
}
\tag{6.4}
$$

## 7. Reproduction and scope

- `searches/internal_singleton_shared_law_contraction.py` computes the exact
  cubic Walsh spectra, theorem coefficient, failed slice diagnostic, ledger,
  next gate, and vertical-mixture stress test.
- `tests/internal_singleton_shared_law_contraction.py` constructs the full
  \(q=4\) cubic endpoint, attains the twirled optimum, tests nonsymmetric
  laws, checks the quintic completion and collision energy, and protects all
  \(q=32\) values.
- `./run_round3_checks.sh` runs the inherited and Round 3 suite.

This theorem does not certify the remaining 872 balanced entries, replace
the known physical lower witnesses by upper bounds, or prove the adaptive
lift.  It does establish the go/no-go orbit below its required gate and
identifies shared-law twirling plus Gram completion as a reusable mechanism
for the next frontier.

## 8. Subsequent fifth-orbit result

The next gate in (6.4) has passed.  For
\((3,1,1,5):(0,0,1,4)\), exact quintic row energy plus rank--Frobenius gives
coefficient \(0.0311889051<1/32\).  With this fifth theorem inserted, the
diagnostic is \(0.3323657941\), leaving \(0.0009675392\); twenty entries are
controlled and 868 remain provisional.  The next reranked orbit is
\((1,1,3,5):(0,1,1,3)\), with gate \(0.0570749885\).  See
`column_cubic_quintic_row_contraction.md`.
