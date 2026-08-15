# Active--passive boundary map

Date: 2026-07-15

## Purpose

The program should not track only a passive lower bound.  This file keeps the
best accepted evidence on both sides of the comparison in one place and
makes missing upper or lower information visible.

## Current rigorous map

| Question | Best accepted result | What it establishes | What it does not establish |
|---|---|---|---|
| Passive asymptotic lower bound | \(D_{\mathsf P}^{\rm hard}=\Omega(N^{1/12})\). | The complete four-initial-state image plus centered repairs close every terminal obstruction from levels six through ten; all sensitive levels five through twelve have decay \(N^{-1}\). | Useful constants at \(N=1024\), an optimizer-excluding physical-frame theorem, or anything close to optimal scaling. |
| Passive finite-size lower bound | No proof yet of passive dose \(>6\) at \(N=1024\).  Ten chain-aware orbit theorems give a coarse completion diagnostic \(0.333132605<1/3\), but 848 balanced entries remain provisional. | The finite-size route remains quantitatively viable, with only \(0.000200728\) diagnostic margin after replacing the tenth provisional charge by its rigorous coefficient. | The diagnostic is not a theorem until all physical coefficients and adaptivity are controlled; the shrinking margin makes a hard-instance or global-contraction pivot increasingly plausible. |
| Passive algorithmic upper evidence | No competitive passive upper-bound frontier has yet been consolidated in Round 3. | This is an explicit missing half of the boundary map. | The absence of an entry is not evidence that the passive lower bound is close to optimal. |
| Active upper bound | \(D_{\mathsf A}^{\rm hard}\le6\). | Six charged traversals suffice for the constant-margin task. | Six is not known to be necessary. |
| Active lower evidence near six | Two complete folded-chain flags have optimal collective error \(0.3611610554>1/3\) on an exact endpoint ensemble at \(N=1024\). | Deleting the third flag or merely decoding two complete flags more cleverly does not reduce the known construction below six. | A general active hard-dose-five lower bound; coherent use of a fifth traversal remains open. |
| Hard-instance robustness | The repaired interpolation witness proves the asymptotic theorem; the signed-permutation witness has promising finite-size structure. | At least two useful witness roles are known. | That either witness is optimal for finite size, constants, or adaptive stability. |
| Mechanism | Reverse-tree dimension decay is proved; the complete terminal image is classified and repaired under assigned, all-projective, and centered expansions. | Passive tensor structure supplies \(N^{-1}\) throughout levels 5--12; a full-mask lower witness proves the limiting level-twelve grouped graph norm is exponent-sharp. | Whether physical posterior frames exclude that optimizer, whether a different hard instance yields more, and which mechanism matters at realistic size. |

## Benchmark questions

The central finite-size benchmark is \(N=1024,D=6\), because it directly
compares with the known active protocol.  When a certificate becomes cheap
enough to evaluate, also report it at nearby powers of two.  This distinguishes
a real finite-size trend from a single optimized checkpoint.

For every benchmark, record four numbers or honest blanks:

| Size | passive lower | passive upper/protocol | active lower | active upper |
|---:|---:|---:|---:|---:|
| asymptotic | \(\Omega(N^{1/12})\) | not yet consolidated | no growing lower bound established here | \(6\) |
| \(N=1024\) | \(>6\) not yet proved | not yet consolidated | two-complete-flag family fails; general \(5\)-dose case open | \(6\) |

Add rows only when the model, error target, promise, and hard-dose meter match
the main problem.  A relaxation or restricted protocol family belongs in the
evidence column of a note, not in a bound cell.

## Route implications

- Closing the signed-permutation one-batch ledger changes only the passive
  finite-size row.  The adaptive lift is still required.
- A stronger reverse-tree contraction changes the asymptotic passive row and
  may also improve finite-size constants if made explicit.
- A passive counterprotocol supplies missing upper evidence and may reveal
  that the current lower-bound target is not close to the true boundary.
- An active five-dose protocol changes the comparison benchmark immediately;
  a restricted-family obstruction does not.
- A better hard instance matters only if it improves the full tester norm,
  conditioning, and adaptive stability under the common scorecard.

## Update rule

Update this map after any theorem, certified protocol, hard-instance
comparison, or obstruction that changes a bound or removes a live route.
For each update, link the supporting artifact and preserve the distinction
between a general result, a restricted family, and a numerical diagnostic.
