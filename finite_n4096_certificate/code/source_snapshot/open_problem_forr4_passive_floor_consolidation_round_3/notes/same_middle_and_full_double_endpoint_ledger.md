# Same-middle contraction and the full double-endpoint ledger

Date: 2026-07-14

Status: the eight same-middle masks of the balanced
\((3,1,1,3)\) obstruction now have exact arbitrary-diagonal coefficients.
The alternating mixed coefficient is also sharpened by an exact vertical
type calculation.  Combining these facts with the complete whole-block cut
table and all 210 dose-six occupation states puts the full 64-cut
double-endpoint ledger below the previously allocated \(N=1024\) margin.
The occupation certificate is conservative floating-point numerics with a
supporting hyperplane; intervalization and the other Fourier profiles remain
open.

## 1. Endpoint orbit data

Put \(N=q^2\), \(n=q-1\), and

$$
A=q^2-2q+2.
\tag{1.1}
$$

For a nonzero pair difference \(x\), let

$$
W_x(i,E)=w_i(E),\qquad \sigma(E)=x,
\tag{1.2}
$$

where \(w_i(E)\) is the endpoint coefficient from the cubic-to-singleton
moment.  Its fixed-singleton squared energy is

$$
S_x=\sum_{\sigma(E)=x}|w_i(E)|^2.
\tag{1.3}
$$

There are two difference classes.  If the hidden-column component of
\(x\) vanishes, call \(x\) vertical; otherwise call it nonvertical.  Direct
orbit counting gives

$$
S_V={A\over2n},\qquad S_\perp={1\over n}.
\tag{1.4}
$$

The endpoint weight matrix has exact nuclear norms

$$
t_V:=\|W_V\|_1=\sqrt2 A,
\qquad
t_\perp:=\|W_\perp\|_1=\sqrt2 q.
\tag{1.5}
$$

One way to verify (1.5) is to diagonalize \(W_x^*W_x\).  In the vertical
case its nonzero singular values are

$$
\sqrt2\quad(A/2\text{ times}),
\qquad
{\sqrt2 A\over2n}\quad(n\text{ times}),
\tag{1.6}
$$

and in the nonvertical case they are

$$
{\sqrt2\over n}\quad(qn/2\text{ times}),
\qquad
\sqrt2\quad(q/2\text{ times}).
\tag{1.7}
$$

The regression reconstructs these spectra directly at \(q=4,8,16\).

## 2. Translation twirl for same-middle masks

Put both middle singleton coordinates on one side of the occurrence cut.
Independent XOR translations in the four physical blocks act by row and
column permutations and diagonal Walsh signs.  Root fidelity is jointly
concave, so twirling cannot lower the weighted nuclear norm.

For equal endpoint orientations, the row law becomes uniform on the two
endpoint singletons.  The column law retains only a joint distribution
\(Q_{x,y}\) on the two pair differences.  Summing the two middle Walsh
coordinates makes the row Gram diagonal and gives exactly

$$
\Phi_{=}(Q)
={2\over N^{3/2}}
\sqrt{\sum_{x,y}Q_{x,y}S_xS_y}.
\tag{2.1}
$$

Concentration on two vertical differences is optimal.  Therefore

$$
\boxed{
\gamma_{=}^{\rm same}
={A\over(q-1)q^3}.
}
\tag{2.2}
$$

For mixed endpoint orientations, the twirled row and column laws retain
separate difference distributions \(P_y,Q_x\).  The column Gram splits by
\(x\), with core \(W_x^*W_x\), and

$$
\Phi_{\ne}(P,Q)
={2\over N^{5/2}}
\sqrt{\sum_yP_yS_y}
\sum_x\sqrt{Q_x}\,t_x.
\tag{2.3}
$$

Maximize the first factor by (1.4) and the second by Cauchy--Schwarz.  Since
there are \(n\) vertical and \(q n\) nonvertical differences,

$$
\boxed{
\gamma_{\ne}^{\rm same}
={2\sqrt{A(A^2+q^3)}\over q^5}.
}
\tag{2.4}
$$

At \(q=32\),

$$
\gamma_{=}^{\rm same}=0.000947029359879\ldots,
\qquad
\gamma_{\ne}^{\rm same}=0.001809666065927\ldots.
\tag{2.5}
$$

Both are far below the formerly charged \(1/q=0.03125\).  Direct matrices
at \(q=2\), with random difference laws, agree with (2.1) and (2.3).

Exactly one endpoint can instead be balanced while the other cubic endpoint
is unsplit.  If the endpoint pair is the only column variable, (1.5) and the
whole-cubic column Bessel energy give

$$
\gamma_{\rm hybrid}^{\rm pair}
={2\sqrt{(q-1)(A^2+q^3)}\over q^5}.
\tag{2.6}
$$

In the transposed same-middle placement, the whole cubic link has weighted
Frobenius norm at most one and rank at most \(N\).  The remaining endpoint
energy is (1.4), so

$$
\gamma_{\rm hybrid}^{\rm whole}
={1\over q^2}\sqrt{{A\over q-1}}.
\tag{2.7}
$$

At \(q=32\), these are \(0.000324856645\) and \(0.005440098115\).
Using the larger value for either orientation is safe and replaces another
previous \(1/q\) charge.  Direct \(q=2\) matrices attain both formulas.

## 3. Sharper alternating mixed coefficient

For the alternating mixed orientation, Section 3 of
`weighted_double_endpoint_contraction.md` reduces the arbitrary-diagonal
problem to diagonally weighted punctured Walsh matrices.  Restrict both
difference laws to the vertical type.  The group
\(GL(\log_2q,2)\) acts by \(x\mapsto Lx\) on one law and
\(y\mapsto L^{-T}y\) on the other.  Averaging any pair of laws under this
joint action makes both uniform, while joint concavity cannot decrease the
norm.  Hence the uniform vertical calculation is the exact arbitrary-law
constant inside this type block.

Let \(r=q/2\), \(h=r-1\), and \(C=q^2-2q+4\).  The needed punctured Walsh
nuclear sums are

$$
T_W=1+(q-2)\sqrt q,
\tag{3.1}
$$

$$
T_C=(r-2)\sqrt q+\sqrt{r+1},
\tag{3.2}
$$

$$
T_0=(q/4-2)\sqrt q+\sqrt{2q+1},
\qquad
T_1=1+(r-2)\sqrt r.
\tag{3.3}
$$

Here \(T_0\) and \(T_1\) correspond to orthogonal and nonorthogonal pairs
of nonzero hyperplane normals.  The exact vertical coefficient is

$$
\gamma_{VV}
={1\over q^6}
\left{
{C^2\over q-1}T_W
+4qCT_C
+4q^2(hT_0+rT_1)
\right}.
\tag{3.4}
$$

At \(q=32\), this is \(0.0202231464018\ldots\).  Replace the old
vertical/vertical entry of the analytic three-type matrix by (3.4), leaving
all eight other entries at their proven rank--Frobenius upper bounds.  An
exact-rational upper for every radical and Collatz--Wielandt with

$$
v=(1,0.054155,0.054752)
\tag{3.5}
$$

give

$$
\boxed{
\gamma_{\ne}^{\rm alt}<0.020343.
}
\tag{3.6}
$$

The earlier bound was \(0.0306880103312\).

## 4. Whole-block cuts

The global occupation optimizer can activate endpoint splits \(0|3\) and
\(3|0\), not only the balanced \(1|2\) masks.  When all four blocks are
unsplit, use the complete three-link cut table rather than charging every
mask \(1/q\).  With endpoint and central coherences all at most \(1/q\),
the safe coefficients are

| whole-block cut, up to complement | coefficient |
|---|---:|
| empty | \(q^{-3}\) |
| any singleton | \(q^{-2}\) |
| adjacent pair touching an endpoint | \(q^{-2}\) |
| central, crossing, or separated pair | \(q^{-1}\) |

The crossing row uses the Bessel-refined path theorem, so no growing outer
operator norm is inserted.

## 5. All 64 cuts and all occupation states

For an occupation distribution \(\rho\) on

$$
\{n\in\mathbb Z_{\ge0}^4:\ |n|\le6\},
\tag{5.1}
$$

define, for a split \(s\le(3,1,1,3)\),

$$
m_s(\rho)
=\mathbb E_\rho
\prod_{b=1}^4\binom{n_b}{s_b}.
\tag{5.2}
$$

If \(\gamma_s\) is the best currently proved coefficient for that cut, the
complete shared-occupation triangle ledger is

$$
F(\rho)
=\sum_s\gamma_s
\sqrt{m_s(\rho)m_{a-s}(\rho)},
\qquad a=(3,1,1,3).
\tag{5.3}
$$

This is concave in \(\rho\).  At any interior candidate \(\rho_0\), its
tangent gives the global upper

$$
F(\rho)\le
\max_{|n|\le6}\nabla F(\rho_0)(n).
\tag{5.4}
$$

CLARABEL supplies a candidate, and direct evaluation of (5.4), with a
conservative numerical allowance, gives

$$
\sup_\rho F(\rho)<0.4987.
\tag{5.5}
$$

After degree-eight attenuation,

$$
\left({5\over6}\right)^8\sup_\rho F(\rho)
<0.1160.
\tag{5.6}
$$

The margin left after the minimal sector and promise conditioning is
\(0.160358131958\).  Thus the formerly dominant double-endpoint profile now
fits by more than

$$
0.0443.
\tag{5.7}
$$

This is a full 64-cut occupation result, unlike the earlier
\((2,1,1,2)\) diagnostic.  It is not a complete passive lower bound: the
other degree-six through degree-twelve profiles must share the same total
TV budget.  The new slack makes their exact joint ledger a plausible next
calculation, but is not itself a bound on those profiles.

Reproduction:

- `searches/same_middle_weighted_bound.py`;
- `searches/mixed_endpoint_weighted_bound.py`;
- `searches/double_endpoint_occupation_optimization.py`;
- `tests/same_middle_weighted_bound.py`; and
- `tests/double_endpoint_occupation_optimization.py`.
