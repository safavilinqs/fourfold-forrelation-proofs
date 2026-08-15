# Alternating double-endpoint spectrum

Date: 2026-07-14

Status: OP3-1 resolved.  The row Gram has an exact two-dimensional Walsh
symbol for every power-of-two \(q\).  Integer orbit counts give a complete
\(q=32\) spectrum certificate.  This is still a uniform-weight result, not
the physical-diagonal theorem required by OP3-2.

## 1. Endpoint coefficient

Write the \(N=q^2\) coordinates as pairs over
\(\mathbb F_2^{\log_2 q}\), and let \(\operatorname{col}(i)\) be the second
coordinate.  For an unordered pair \(E=\{u,v\}\), put
\(\sigma(E)=u\oplus v\).  Direct signed-permutation averaging gives

$$
A(i;E,b)=w_i(E)H_N(i\oplus\sigma(E),b),
\tag{1.1}
$$

where \(w_i(E)=0\) if \(i\in E\), and otherwise

$$
w_i(E)=
\begin{cases}
1,
&|\{\operatorname{col}(i),\operatorname{col}(u),
       \operatorname{col}(v)\}|=1,\\[2mm]
-1/(q-1),
&|\{\operatorname{col}(i),\operatorname{col}(u),
       \operatorname{col}(v)\}|=2,\\[2mm]
0,
&|\{\operatorname{col}(i),\operatorname{col}(u),
       \operatorname{col}(v)\}|=3.
\end{cases}
\tag{1.2}
$$

This formula reproduces the independently enumerated endpoint tensor at
\(q=2,4\).

## 2. Exact block reduction

The alternating lift is

$$
K((i,b,d),(E,c,F))
=A(i;E,b)H_N(b,c)A(d;F,c).
\tag{2.1}
$$

Apply the normalized Walsh transform to \(b\), with dual label \(t\).
Walsh orthogonality gives

$$
\widetilde K((i,t,d),(E,c,F))
={w_i(E)w_d(F)\over N}\,
 \mathbf 1\{c=t\oplus i\oplus\sigma(E)\}
 \chi_{d\oplus\sigma(F)}(c).
\tag{2.2}
$$

Two rows can therefore have nonzero inner product only when

$$
t\oplus i=t'\oplus j.
\tag{2.3}
$$

The row Gram splits into \(N\) blocks indexed by this conserved label.  The
blocks are diagonally unitarily equivalent.  After removing the harmless
phase, their common core is

$$
T((i,d),(j,\ell))
={1\over N^2}S(d\oplus\ell)
 B(i\oplus j,d\oplus\ell),
\tag{2.4}
$$

where

$$
S(r)=\sum_Ew_0(E)w_r(E),
\qquad
B(h,r)=\sum_Ew_0(E)w_h(E)\chi_{\sigma(E)}(r).
\tag{2.5}
$$

Thus \(T\) is an XOR convolution on two copies of the \(N\)-element group.
Its full spectrum is its two-dimensional Walsh transform.

## 3. Integer certificate

Let \(U_i(E)=(q-1)w_i(E)\in\{q-1,-1,0\}\), and define the integer table

$$
C(h,x)=
\sum_{\substack{E\ {\rm unordered}\\\sigma(E)=x}}
U_0(E)U_h(E).
\tag{3.1}
$$

Put

$$
s(r)=\sum_x C(r,x),
\qquad
b(h,r)=\sum_xC(h,x)\chi_x(r).
\tag{3.2}
$$

For Walsh frequencies \((\alpha,\beta)\), the block eigenvalues are exactly

$$
\boxed{
\lambda_{\alpha,\beta}
={1\over N^2(q-1)^4}
\sum_{r,h}s(r)b(h,r)\chi_\alpha(r)\chi_\beta(h).
}
\tag{3.3}
$$

Every operation in (3.1)--(3.3) is an integer count or integer Walsh
transform before the displayed denominator is applied.  The uniform
normalized nuclear coefficient is

$$
\gamma_q
={1\over N\binom N2}
\sum_{\alpha,\beta}\sqrt{\lambda_{\alpha,\beta}}.
\tag{3.4}
$$

The implementation checks positivity, full multiplicity \(N^2\), and the
exact block trace

$$
\sum_{\alpha,\beta}\lambda_{\alpha,\beta}
=\left({N+2\over2}\right)^2.
\tag{3.5}
$$

## 4. Results

| \(q\) | \(N\) | spectral classes | block rank | \(\gamma_q\) |
|---:|---:|---:|---:|---:|
| 2 | 4 | 3 | 16 | 0.4715918158911433 |
| 4 | 16 | 19 | 256 | 0.06420087162467489 |
| 8 | 64 | 20 | 4096 | 0.01117494129649916 |
| 16 | 256 | 20 | 65536 | 0.002004962418830680 |
| 32 | 1024 | 20 | 1048576 | 0.0003575935171398254 |

The \(q=2,4\) values match the independent direct row-Gram calculation.
At \(q=32\), the common denominator in (3.3) is

$$
968381956096.
$$

The twenty exact numerator/multiplicity pairs are printed by the search
artifact.  As compact independent checks, the smallest numerator is
\(889351684\) with multiplicity \(15872\), the largest is
\(7945127554564\) with multiplicity \(15376\), and the weighted numerator
sum equals the denominator times \(513^2=263169\).

## 5. What failed in the first reduction

Writing a pair as \(\{u,u\oplus x\}\) is two-to-one:

$$
\{u,u\oplus x\}=\{u\oplus x,u\}.
$$

Equivalently, translation by \(x\) stabilizes the unordered pair.  The pair
orbit is a quotient by \(\{0,x\}\), not a free \(N\)-element translation
orbit.  The first candidate autocorrelation reduction did not retain this
stabilizer consistently, which is why it lost exactly \(N\) directions per
block:

| \(q\) | old candidate coefficient | exact coefficient | old block rank | exact block rank |
|---:|---:|---:|---:|---:|
| 2 | 0.426776695297 | 0.471591815891 | 12 | 16 |
| 4 | 0.063628751110 | 0.064200871625 | 240 | 256 |

Equation (3.1) avoids a choice of representatives and counts every unordered
pair once.  It is therefore insensitive to this bookkeeping trap.

## 6. Consequence for the finite-size route

This matrix has the same \(1|2\) orientation at both decorated endpoints.
Its old local row/column coefficient at \(q=32\) is
\(E_2=0.03030494\), not the mixed-orientation worst value
\(0.12321552\).  The uniform coefficient is much smaller, but nonuniform
weights can partially undo that cancellation.

The translation-twirled arbitrary-diagonal optimization is now below
\(0.010905\).  The distinct mixed \(1|2\) endpoint orientation has the
rationally certified analytic bound \(0.020343\).  See
weighted_double_endpoint_contraction.md and
same_middle_and_full_double_endpoint_ledger.md for the corrected comparison
and complete occupation ledger.

Reproduction:

- searches/alternating_double_endpoint_spectrum.py;
- tests/alternating_double_endpoint_spectrum.py; and
- the independent direct \(q=2,4\) benchmark in the frozen round-two folder.
