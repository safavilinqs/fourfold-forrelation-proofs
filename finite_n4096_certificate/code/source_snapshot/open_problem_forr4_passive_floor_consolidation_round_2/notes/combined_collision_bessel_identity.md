# Combined collision and insertion/Bessel identity

Date: 2026-07-14

Status: explicit local bookkeeping identity; this closes the question
whether the parity-support collision factor duplicates the insertion mass.

## 1. One amplitude entry

Let \(q\) be a probability distribution on parity supports \(S\), with
\(|S|\le t\). Put \(a\) labeled marks in this amplitude entry. For an
ordered tuple \(x=(x_1,\ldots,x_a)\) of distinct marked coordinates and a
base support \(R\) disjoint from \(x\), define

$$
\tau(R,x)=R\cup\{x_1,\ldots,x_a\}.
$$

For a fixed base \(R\) and output support \(S\), the preimage consists of
the orderings of \(S\setminus R\). Hence

$$
\max_{R,S}|\{x:\tau(R,x)=S\}|\le a!.
\tag{1.1}
$$

This is the collision multiplicity used by the collision-aware frame
lemma. It does not sum over insertion choices.

Now sum the diagonal frame mass over all bases and ordered marked tuples:

$$
\begin{aligned}
\sum_{\substack{R,x\\R\cap\{x_i\}=\varnothing}}
q(\tau(R,x))
&=\sum_S q(S)(|S|)_a\\
&\le t^a.
\end{aligned}
\tag{1.2}
$$

Every pair \((R,x)\) on the left is uniquely obtained from \(S\) and an
ordered choice of \(a\) distinct elements of \(S\). Thus (1.2) is an
identity before the dose inequality is applied.

## 2. Bilateral entry and graph fibers

With \(a\) ket marks and \(b\) bra marks, the two complete-frame square
sums are bounded by (1.2) with powers \(t^a\) and \(t^b\). The
collision-aware bilateral lemma therefore costs

$$
\sqrt{a!b!}\,t^{(a+b)/2}.
\tag{2.1}
$$

Coordinates assigned to opposite-entry graph fibers are outer indices of
the sliced covariance family. Summing them replaces its constant diagonal
by

$$
\sum_s\operatorname{diag}(M_sM_s^*)
=N^{|B|-e}\mathbf1.
$$

This sum multiplies (1.2); it does not create another choice of \(R\) or
\(x\). Local Walsh character factors from the distinct-label repair are
diagonal unitaries and leave the covariance diagonal unchanged.

Therefore collision multiplicity, graph-fiber cancellation, and insertion
mass occur in one local square sum. There is one falling-factorial charge,
not two.

## 3. Scope

Repeated physical occupations have already been quotiented to the odd
parity support \(S\). Equation (1.2) depends only on \(|S|\le t\), so
repeated modes and vacuum support require no separate count.

The exact finite enumeration in
tests/insertion_collision_ledger.py checks (1.1)--(1.2) for every support in
small universes and random rational support laws.
