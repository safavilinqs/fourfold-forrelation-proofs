# Weighted endpoint contraction for the quadratic-bent plant

Date: 2026-07-14

Status: proved for every endpoint link \(M_{a,1}\) and its transpose,
uniformly in support degree and dimension.

## 1. Exact XOR-row form

Let \(X=(-1)^{q+\ell}\) be sampled uniformly from the
nondegenerate-quadratic orbit, including its uniform linear shift
\(\ell\), and let \(Y=H_NX\).  For a distinct-coordinate set
\(A\subseteq\mathbb F_2^n\), define

$$
M_{a,1}(A,y)=\mathbb E[X_AY_y],
\qquad a=|A|.
$$

Expand one Walsh coordinate:

$$
Y_y=\sum_zH_{y,z}X_z.
$$

The average over the linear shift vanishes unless

$$
\bigoplus(A\mathbin\triangle\{z\})=0,
$$

which forces \(z=\bigoplus A\).  Therefore

$$
\boxed{
M_{a,1}(A,y)=\mu_AH_{\oplus A,y},
}
\tag{1.1}
$$

where

$$
\mu_A=\mathbb E[X_{A\mathbin\triangle\{\oplus A\}}],
\qquad |\mu_A|\le1.
$$

The potentially large family of sets sharing one XOR produces repeated
Walsh rows, but no entry larger than \(N^{-1/2}\).

## 2. Physical diagonal weighting

Let \(p_A,q_y\ge0\), with total masses

$$
P=\sum_Ap_A,\qquad Q=\sum_yq_y.
$$

Equation (1.1) gives

$$
\begin{aligned}
\|D_p^{1/2}M_{a,1}D_q^{1/2}\|_F^2
&=\sum_{A,y}p_Aq_y|\mu_A|^2|H_{\oplus A,y}|^2\\
&\le {PQ\over N}.
\end{aligned}
$$

The matrix has at most \(N\) columns, hence rank at most \(N\).  The
rank--Frobenius inequality proves

$$
\boxed{
\|D_p^{1/2}M_{a,1}D_q^{1/2}\|_1
\le\sqrt{PQ}.
}
\tag{2.1}
$$

The same result holds for \(M_{1,a}\) by transposition and invariance of
the orbit under Walsh transform.

For \(a=5\), the sharper multiplier bound
\(|\mu_A|\le2/(N-2)\) gives the \(O(1/N)\) estimate recorded in
notes/quadratic_bent_candidate.md.  For \(a=7,9,\ldots\), some
\(\mu_A\) equal one because affine-flat parity relations are deterministic;
(2.1) remains dimension-free and sharp enough to show that these
relations do not amplify the endpoint link beyond the minimal \(M_{1,1}\)
scale.

## 3. Scope

This closes every outer link having a singleton on one side.  It does not
bound \(M_{a,b}\) when \(a,b\ge2\), especially the internal link in
profiles such as \((1,5,5,1)\).  Those weighted two-sided sectors are now
the smallest unresolved association-scheme calculation for the
quadratic-bent plant.

Reproduction: tests/quadratic_endpoint_weighted_bound.py.
