# Adjacent double-cubic occurrence contraction

Date: 2026-07-14

Status: the record-one sector of \((3,3,1,1)\), and its path reversal,
is bounded for every fixed unordered occurrence split and arbitrary
correlated diagonal weights.  A uniform incidence bound gives coefficient
at most \(0.009277\) at \(q=32\), or \(0.002158\) after degree-eight
attenuation.  The joint sum over occurrence assignments is not included.

## 1. The endpoint-to-L link

The endpoint record-one cubic supports have two forms:

1. three cells in one column;
2. two cells in one column and one cell in another column.

The compatible middle cubic supports are the three-cell L-shapes in a
two-by-two rectangle.  For an endpoint support \(A\) and L-shape \(B\),
the squared signed-permutation link entry has baseline value

$$
a_0={1\over q^2(q-1)^2}
$$

and, when the two Walsh orthogonality conditions hold, exceptional value

$$
a_1={ (q+2)^2\over q^2(q-1)^2(q-2)^2}.
\tag{1.1}
$$

These are the only two values.  The exceptional condition is an ordinary
incidence relation between the row difference of \(B\), the column
difference of \(A\), and the two endpoint Walsh labels; it introduces no
extra labeled-mark multiplicity.

## 2. Uniform slice bound

Let \(D^E_k\) be the maximum number of endpoint supports containing a
fixed \(k\)-cell subset, and let \(D^L_\ell\) be the analogous L-shape
degree.  Direct incidence counting gives

$$
\begin{aligned}
D^E_0&=q\binom q3+q^2(q-1)\binom q2,\\
D^E_1&=\binom{q-1}2+(q-1)\binom q2+q(q-1)^2,\\
D^E_2&=q^2-2,\qquad D^E_3=1,
\end{aligned}
\tag{2.1}
$$

and

$$
D^L_0=q^2(q-1)^2,\quad
D^L_1=3(q-1)^2,\quad
D^L_2=2(q-1),\quad D^L_3=1.
\tag{2.2}
$$

Therefore the squared \(M_{3,3}\) slice obtained by fixing \(k\) endpoint
and \(\ell\) middle occurrences is at most

$$
S_{k\ell}\le a_1D^E_kD^L_\ell.
\tag{2.3}
$$

For the path \(M_{3,3}M_{3,1}H_N\), summing the singleton neighbor of the
L-shape contributes \(1/(q-1)^2\) to the row energy.  Reversing the
flattening and using the pointwise singleton-link bound contributes
\(1/[q^2(q-1)^2]\) to the column energy.  We deliberately discard any
additional gain from the final Hadamard link.  Thus the split
\((k,\ell)\) has weighted coefficient at most

$$
\boxed{
\Gamma_{k\ell}(q)
\le
\min\left\{
{\sqrt{a_1D^E_kD^L_\ell}\over q-1},
{\sqrt{a_1D^E_{3-k}D^L_{3-\ell}}\over q(q-1)}
\right\}.
}
\tag{2.4}
$$

This follows from the two row/column Schatten factorizations and allows
arbitrary correlated diagonal masses within the fixed split.  The
singleton placements that make a link internal only improve the bound;
the unsplit endpoint cases are also covered by the earlier whole-block
Bessel contraction.

At \(q=32\), evaluating all sixteen elementary expressions in (2.4)
shows that the largest is the \((k,\ell)=(2,2)\) split:

$$
\boxed{
\max_{k,\ell}\Gamma_{k\ell}(32)
\le 0.009276958.
}
\tag{2.5}
$$

Multiplication by \((5/6)^8\) gives \(0.002157524\).

## 3. Exact check of the worst slice

The uniform estimate uses the exceptional value for every compatible
pair.  Counting the exceptional pairs exactly for the \(2|1,2|1\) split
gives

$$
F_{22}={2(q^2-2)\over q^2(q-1)}
       +{8\over(q-1)(q-2)},
\tag{3.1}
$$

while the complementary \(1|2,1|2\) slice has

$$
G_{11}=9+
{3(q-1)\{q-2+3q(q-1)\}\over2q^2}.
\tag{3.2}
$$

Consequently this particular split has coefficient

$$
\min\left\{{\sqrt{F_{22}}\over q-1},
{\sqrt{G_{11}}\over q(q-1)}\right\},
\tag{3.3}
$$

which equals \(0.00871519\) at \(q=32\).  Complete orbit enumeration at
\(q=4\) and \(q=8\) verifies (3.1)--(3.2) and independently finds this
split to be the largest of all sixteen.  The rigorous \(q=32\) statement,
however, uses the uniform table (2.4), not extrapolation from those finite
checks.

## 4. Scope

This closes the local fixed-split occurrence problem for the record-one
\((3,3,1,1)\) sector and, by reversal, \((1,1,3,3)\).  Together with the
double-endpoint calculation, every formerly hard record-one degree-eight
profile now has a small fixed-split coefficient.

The unresolved step is common to all of them: the different split masks
must be packed jointly against one passive occupation law.  No triangle
sum over occurrence assignments has been charged here.

Reproduction: `searches/adjacent_double_cubic_orbit_scaling.py` and
`searches/adjacent_double_cubic_slice_energies.py`.
