# Realistic-size program: exact signed-permutation plant

Date: 2026-07-14

Target: prove that every passive protocol of hard dose at most six has transcript total variation below $1/3$ at $N=1024$.

## 1. Why the Gaussian asymptotic route is not viable

The repaired bound

$$
\operatorname{TV}\le C(1+D)^{12}/\sqrt N
$$

is correct asymptotically but gives $C\,7^{12}/32\approx4.33\times10^8C$ at $D=6$, $N=1024$. Constant extraction or expository tightening cannot bridge this gap.

## 2. Candidate hard instance

Use the exact Maiorana--McFarland signed-permutation plant from Note 14. For $N=q^2$ with $q$ a power of two, it gives distributions $Q_+$ and $Q_-$ supported pointwise on

$$
F_{4,H}=+1\quad\text{and}\quad F_{4,H}=-1.
$$

At $N=1024$, $q=32$. There is no conditioning loss and no Gaussian interpolation tail. Hypothesis sensitivity requires a matching event through each of three independent hidden signed permutations.

## 3. Exact finite occurrence budget

For two blocks containing at most $a$ and $b$ ket/bra occurrences, define

$$
L_q(a,b)=
\sum_{r=1}^{\min(a,b)}
{\binom ar\binom br\over\binom qr}.
$$

This is the union weight obtained by choosing an odd label subset of size $r$ on each side and asking a uniform hidden permutation to map one set to the other.

If a paired amplitude history has block occurrence counts $(n_1,n_2,n_3,n_4)$, the three independent match weights multiply:

$$
B_q(n_1,n_2,n_3,n_4)
=L_q(n_1,n_2)L_q(n_2,n_3)L_q(n_3,n_4).
$$

A hard-dose-$D$ ket/bra pair has

$$
n_1+n_2+n_3+n_4\le2D.
$$

Exact enumeration for $q=32$, $D=6$ gives

$$
\max B_{32}
=B_{32}(2,4,4,2)
={17497415\over442336768}
\approx0.0395568.
$$

Including the factor two between the positive and negative hypotheses gives

$$
2\max B_{32}\approx0.0791135<1/3.
$$

Thus a recording/composition theorem may lose a further factor as large as about $4.21$ and still rule out dose six at $N=1024$.

The same budget is not sufficient at $N=256$ ($q=16$), where it is about $0.804$ for dose six. At $N=4096$ it is about $0.00879$. The $N=1024$ target sits in a useful, nontrivial range.

## 4. Exact theorem to prove

The preferred target is a complete-outcome recording inequality of the form

$$
\operatorname{TV}(Q_+^{\mathcal T},Q_-^{\mathcal T})
\le K\,
2\max_{n_1+\cdots+n_4\le2D}
B_q(n_1,n_2,n_3,n_4),
$$

for every classically adaptive passive tree, with a proved constant $K<4.21$; ideally $K=1$.

The statement must be proved at the level of the physical four sign blocks. It may not reveal the hidden permutations or grant oracle access to their factors.

## 5. Proof route

1. Expand a signed terminal transcript functional into paired ket/bra occupation histories. The total number of physical phase occurrences is at most $2D$ on every branch pair.
2. For each of the three hidden signed permutations, expose only the odd label sets transmitted between its $X$ and $Y$ appearances.
3. Apply the exact moment formula: a fixed size-$r$ record vanishes unless the two sets have equal size, and then has magnitude at most $1/\binom qr$.
4. Sum record choices within each block to obtain $L_q$.
5. Use independence of the three permutations to multiply the link weights.
6. Prove that reverse POVM completeness turns the remaining adaptive amplitude coefficients into a positive subprobability or a norm-one recording operator, with loss $K<4.21$ and no dependence on outcome width or depth.

Step 6 is the real problem. Entrywise absolute summation is too weak: the averaged matching matrix can have operator norm one even when each entry is $1/\binom qr$. The proof must use the product-mask structure and passive complete frames, not only the combinatorial event probability.

## 6. Falsification program

- Build the exact fixed-versus-fresh signed-permutation transcript symbol for $q=2,4$ and low dose.
- Use seesaw optimization to search passive root/child frames that exceed the proposed recording budget.
- Formulate constrained tester SDPs for fixed batch structures; ordinary PPT or separability without preparation/POVM normalization is not sufficient.
- Test the dangerous adjacent-pair and singleton--triple placements first, since generic operator norms lose the matching suppression there.
- If a counterexample exploits coherent superposition over many matching records, refine the target to a three-link noncommutative square function rather than an entrywise union bound.

## 7. Backup contraction route

If the exact recording inequality fails, retain the full component factors in the repaired global dichotomy and seek a temporal square-function bound across marked times. To be useful at $N=1024$, the resulting explicit degree-four contribution must be below $1/3$ at $D=6$; an unspecified $D^{O(1)}$ bound is no longer an acceptable milestone.
