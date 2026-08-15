# Gram-dressed tail contraction

Date: 2026-07-14

Status: proved for arbitrary diagonal weights.  This closes every
one-batch support profile \((a,1,1,1)\) and \((1,1,1,a)\), including all
internal ket/bra splits of the decorated endpoint support, at the sharp
coefficient \(N^{-1/2}\).

## 1. Dressing lemma

Let \(T(u,v)\) be a base kernel satisfying

$$
\|D_P^{1/2}TD_Q^{1/2}\|_1
\le\gamma\sqrt{PQ}
\tag{1.1}
$$

for every nonnegative diagonal weight.  Repeat its rows and columns by
arbitrary ancillary indices, and let

$$
G(i,j)=\mathbb E_\omega
[\overline{f_i(\omega)}g_j(\omega)]
\tag{1.2}
$$

be any cross Gram kernel of unit-modulus features on the enlarged row and
column indices.  Then

$$
\boxed{
\|D_p^{1/2}(G\circ\widetilde T)D_q^{1/2}\|_1
\le\gamma\sqrt{pq}.
}
\tag{1.3}
$$

Indeed,

$$
G\circ\widetilde T
=\mathbb E_\omega[D_{f(\omega)}^*\widetilde T D_{g(\omega)}].
$$

The diagonal phases commute with the physical weights and do not change
the trace norm.  For fixed \(\omega\), repeated rows and columns of
\(\widetilde T\) compress isometrically after their diagonal masses are
summed.  Equation (1.1) applies to those aggregated masses.  Convexity of
the trace norm completes the proof.

This lemma permits arbitrary dependence of \(f_i\) and \(g_j\) on every
index on their respective physical sides.  No product assumption on the
probe weight is made.

## 2. The singleton two-link tail

For the last two minimal Hadamard links, consider any ket/bra assignment
of singleton coordinates \(b,c,d\).

If \(b,c\) lie on the same side, \(H_N(b,c)\) is an internal diagonal
factor of magnitude \(N^{-1/2}\); the remaining cross link is a Gram
contraction.  The same argument applies if \(c,d\) lie on the same side.

The only other assignment has \(b,d\) on one side and \(c\) on the other.
Its base kernel is the two-link wedge

$$
T((b,d),c)=H_N(b,c)H_N(c,d).
\tag{2.1}
$$

It has rank at most \(N\) and constant entry magnitude \(1/N\).  Weighted
rank--Frobenius therefore gives

$$
\|D_P^{1/2}TD_Q^{1/2}\|_1
\le {1\over\sqrt N}\sqrt{PQ}.
\tag{2.2}
$$

Uniform weights attain (2.2), so the coefficient is sharp.

## 3. Arbitrarily decorated endpoint

Let the first block have any odd support degree \(a\), while the other
three blocks are singletons.  Split the first-block support arbitrarily
between ket and bra, and place \(b\) on its chosen physical side.  The
first signed-permutation link contributes

$$
G=\mathbb E[(X_{S_1}Y_{S_2})
             (X_{T_1}Y_{T_2})],
\tag{3.1}
$$

which is exactly a unit-feature cross Gram kernel on the two enlarged
physical index families.  The rest of the path is one of the tail kernels
in Section 2.  Applying (1.3) proves

$$
\boxed{
\text{every occurrence split of }(a,1,1,1)
\text{ has weighted coefficient at most }{1\over\sqrt N}.
}
\tag{3.2}
$$

Path reversal proves the same statement for \((1,1,1,a)\).  The cubic
whole-block orbit in
notes/signed_permutation_decoration_compatibility.md attains (3.2), so a
better dimension power is false.

## 4. Remaining boundary

The lemma does not close a decoration in block two or three, because that
physical support participates in two planted links and cannot be absorbed
into one cross Gram dressing of an already-contracting singleton tail.
The next target is a two-sided dressing lemma for profiles
\((1,a,1,1)\), beginning with the exact cubic L-shape sectors.

Reproduction: tests/gram_dressed_tail_contraction.py.
