# Exact one-batch certificate for the minimal chain

Date: 2026-07-14

Status: proved for the minimal four-vertex moment tensor, one passive batch, hard dose at most six, and $N=1024$. This is not yet a certificate for adaptive protocols or higher signed-permutation sectors.

## 1. Exact frame masses

Let $q$ be the diagonal probe distribution and let

$$
n_b(S)=|S\cap\text{block }b|,
\qquad
n_1(S)+\cdots+n_4(S)\le6.
$$

For a ket/bra split $A\subseteq\{1,2,3,4\}$ of the four distinct minimal-chain marks, complete-frame packing gives the exact square masses

$$
m_A=\mathbb E_q\prod_{b\in A}n_b,
\qquad
m_{A^c}=\mathbb E_q\prod_{b\notin A}n_b.
$$

The split contributes at most $\sqrt{m_Am_{A^c}}$ times its graph cut norm.

## 2. Exact cut weights

Let $r_A$ be the binary cross-adjacency rank of the three-edge chain across $A\mid A^c$. The weighted nuclear estimate is

$$
N^{(r_A-3)/2}.
$$

Factoring out the worst scale $N^{-1/2}$ leaves weight

$$
w_A=N^{(r_A-2)/2}.
$$

The four oriented crossing masks have $w_A=1$. The ten nonempty/nonfull singleton, triple, and adjacent masks have $w_A=N^{-1/2}=1/32$. The empty and full masks have $w_A=N^{-1}=1/1024$.

Thus the exact occupation optimization is

$$
F(q)=\sum_{A\subseteq[4]}w_A\sqrt{m_Am_{A^c}}.
$$

This is concave and positively homogeneous in the occupation distribution $q$.

## 3. Optimizer and exact value

Put equal mass $1/4$ on

$$
(1,2,2,1),\quad
(2,1,1,2),\quad
(1,2,1,2),\quad
(2,1,2,1).
$$

Direct evaluation gives

$$
\boxed{
F_6={2337\over256}+{3\sqrt2\over8}
\approx9.65923633589.}
$$

To certify global optimality, differentiate $F$ at this distribution. For every one of the 210 integer occupation states $n\in\mathbb Z_{\ge0}^4$ with $\sum_bn_b\le6$, the directional gradient is at most $F_6$. Equality holds exactly at the four displayed states. The next-largest gradient is lower by

$$
{33\over512}+{3\sqrt2\over32}>0.
$$

Concavity and homogeneity now give, for every distribution $q$,

$$
F(q)\le F(q_*)+\nabla F(q_*)\cdot(q-q_*)
=\nabla F(q_*)\cdot q
\le F_6.
$$

The exact enumeration and sign checks are reproduced by tests/one_batch_minimal_chain_certificate.py.

## 4. Finite-size consequence

The total minimal-chain contribution of any single passive dose-six batch is at most

$$
{F_6\over\sqrt{1024}}
={F_6\over32}
\approx0.301851135497
< {1\over3}.
$$

The margin is

$$
{1\over3}-{F_6\over32}
\approx0.0314822.
$$

This result shows that exact block occupations and cut ranks can reach the desired realistic-size scale; the coarse $D^4/\sqrt N$ ledger cannot.

## 5. Remaining scope

The margin is only about $0.0315$. To prove passive hard dose greater than six, the following still have to be controlled jointly rather than added with coarse triangle inequalities:

1. multiple adaptive batches whose doses sum to six;
2. all higher signed-permutation moment sectors with total Fourier support at most twelve; and
3. sector mixing under outcome-selected complete frames.

The one-batch minimal-chain certificate is therefore a necessary benchmark and a sharp constant calculation, not the final lower bound.
