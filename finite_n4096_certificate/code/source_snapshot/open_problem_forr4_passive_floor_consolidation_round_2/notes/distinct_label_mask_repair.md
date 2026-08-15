# Distinct-label mask repair

Date: 2026-07-14

Status: exact repair of the omitted cross-component Fourier-support mask.

## 1. The omitted mask

The interpolation vertices are distinct marked coordinates. Therefore, if
the vertices in layer (r) are (V_r), the graph summation contains

$$
\Delta_r(x)=\prod_{\{u,v\}\subseteq V_r}{\bf1}\{x_u\ne x_v\}.
$$

Without these factors a disconnected graph tensor is a vertical product of
its component tensors. With them it need not be. This matters even when
every component is singleton in every physical entry.

For two four-layer chains at (N=4), place the second chain in physical
entries ((0,1,3,2)) relative to the first. The exact product-vector witness
in `tests/all_singleton_masked_graph_norm.py` has contraction

$$
{15\over32\sqrt3}>{1\over4}.
$$

The unmasked vertical-product value is (1/N=1/4). Thus literal
multiplicativity after imposing distinct labels is false.

## 2. Character expansion

Identify the Sylvester labels with the group

$$
\Gamma=\mathbb F_2^d,
\qquad |\Gamma|=N.
$$

Character orthogonality gives

$$
{\bf1}\{x=y\}={1\over N}\sum_{\chi\in\widehat\Gamma}
\chi(x)\overline{\chi(y)}.
$$

Consequently

$$
{\bf1}\{x\ne y\}
=\left(1-{1\over N}\right)
-{1\over N}\sum_{\chi\ne1}\chi(x)\overline{\chi(y)}.
\tag{2.1}
$$

The coefficient \(\ell_1\)-mass of (2.1) is

$$
2\left(1-{1\over N}\right)\le2.
$$

If

$$
P_G=\sum_r {|V_r|\choose2},
$$

expanding every factor in every \(\Delta_r\) gives a sum with total
coefficient mass at most \(2^{P_G}\). Since (v(G)\le12), this is an
absolute diagram constant.

## 3. Why every expanded term has the old norm

Each summand of the expansion is a product of one character factor at each
incident vertex. Absorb all factors incident to a vertex into a diagonal
unitary on that vertex's copy of \(\mathbb C^N\). Thus the summand is the
vertical product of the original component tensors after local diagonal
unitaries have been applied.

Local unitaries preserve:

1. every component injective norm and natural-cut operator norm;
2. every sliced covariance diagonal used in the all-assigned case; and
3. every Hilbert preimage norm in the reverse induction.

The Cochrane/Derksen vertical-product theorem therefore applies term by
term. Triangle inequality is taken only over the mask expansion and costs
at most (2^{P_G}), which is absorbed into (C_G). In particular, in the
all-singleton case,

$$
\|\Delta_G\,\mathop{\boxtimes}_C T_C\|_\varepsilon
\le 2^{P_G}\prod_C\|T_C\|_\varepsilon
\le 2^{P_G}N^{-s_0/2}.
\tag{3.1}
$$

The same preliminary expansion makes the all-assigned sliced maps genuine
tensor products again, with unchanged diagonal estimates. Hence the mask
does not alter the (N^{-1/2}) theorem, but the diagram constant is
essential and exact multiplicativity of the masked tensor must not be
claimed.

## 4. Audit consequence

Every use of vertical component multiplicativity or tensorized sliced
covariances must occur after (2.1) is expanded. Treating the distinct-label
mask as though it preserved a vertical product is invalid; the exact
(N=4) witness detects that mistake.
