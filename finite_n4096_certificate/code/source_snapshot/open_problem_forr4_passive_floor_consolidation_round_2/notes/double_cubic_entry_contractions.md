# Double-cubic entry contractions

Date: 2026-07-14

Status: three of the six degree-eight double-cubic profiles have explicit
fixed-occurrence diagonal bounds.  Two separated endpoint/middle profiles
have coefficient \(q/(q-1)^2\).  The adjacent-middle profile has an
additional record-three sector bounded by \(q^2/\binom q3\).  Joint
occurrence packing is not yet done.  This is a historical partial-result
note: the other three fixed-split profiles were subsequently bounded in
the endpoint-occurrence notes.

## 1. Separated endpoint and middle decorations

Consider \((3,1,3,1)\).  The third-block cubic support lies between two
singleton records, so it is an L-shape.  Its two link moments each have
magnitude at most \(1/[q(q-1)]\).  The first endpoint-to-singleton link
has magnitude at most \(1/q\).  Hence

$$
|K|\le {1\over q^3(q-1)^2}.
\tag{1.1}
$$

There are eight marks.  Rank--Frobenius therefore gives, for every fixed
occurrence split,

$$
\boxed{
\|D_p^{1/2}KD_q^{1/2}\|_1
\le {q\over(q-1)^2}\sqrt{PQ}.
}
\tag{1.2}
$$

Path reversal gives the same result for \((1,3,1,3)\).  At \(q=32\),
the coefficient is \(32/961\approx0.033299\).

## 2. Two adjacent middle decorations

For \((1,3,3,1)\), the central signed-permutation record can have size one
or three.

If it has size one, both cubic supports are L-shapes and (1.1)--(1.2)
apply unchanged.

If it has size three, the two outer singleton links can each have full
coherence \(1/q\), while the central record costs
\(1/\binom q3\).  Thus

$$
|K_{r=3}|\le {1\over q^2\binom q3}.
\tag{2.1}
$$

The same rank bound gives

$$
\boxed{
\|D_p^{1/2}K_{r=3}D_q^{1/2}\|_1
\le {q^2\over\binom q3}\sqrt{PQ}.
}
\tag{2.2}
$$

At \(q=32\), this coefficient is

$$
{1024\over4960}\approx0.206452.
$$

This is finite and decays as \(6/q\), but it is too large to discard in a
dose-six triangle sum.  Its pure record-three compound norm should improve
the constant beyond the entry-rank estimate.

## 3. Degree-eight core that remained at this checkpoint

The three profiles not covered here are

$$
(3,3,1,1),\quad(1,1,3,3),\quad(3,1,1,3).
$$

The first pair couples a decorated endpoint directly to a decorated middle
block.  The last decorates both endpoints around a singleton middle edge.
Their dangerous record-one sectors lack the two explicit
\(1/(q-1)\) factors used above and require a Bessel/compound contraction,
not entrywise rank alone.

Reproduction: searches/signed_permutation_double_cubic_entries.py.
