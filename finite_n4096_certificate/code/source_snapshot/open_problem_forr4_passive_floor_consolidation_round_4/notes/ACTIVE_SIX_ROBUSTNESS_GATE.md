# Active six-dose robustness gate

Date: 2026-07-16

Status: exact scalar tolerance for the three-flag majority protocol. This is not a device-specific loss/noise theorem and does not establish experimental feasibility.

## Exact contrast threshold

Let $\mu$ be the signed expectation of each independent binary flag at the promise boundary, after assigning a fair random sign to every no-click or erased event. Its probability of being correct is

$$
p=\frac{1+\mu}{2}.
$$

The three-flag majority error is

$$
\varepsilon_3(\mu)
=\left(\frac{1-\mu}{2}\right)^2(2+\mu).
$$

Solving $\varepsilon_3(\mu)=1/3$ gives

$$
\mu_*=0.226073713789.
$$

The ideal promise-boundary expectation is $1/4$, so the active protocol needs a total multiplicative contrast strictly above

$$
g_*=4\mu_*=0.904294855157.
$$

At unit contrast, the corresponding worst-signed additive expectation-error budget is only

$$
\frac14-\mu_*=0.0239262862108.
$$

## Physical interpretation

For a symmetric scalar model, let

$$
g=\xi\sqrt{T_L T_R}\,\nu\cos\delta,
$$

where $\xi$ is final detection/retention efficiency, $T_L,T_R$ are the two arm power transmissions, $\nu$ is task-mode overlap, and $\delta$ is differential phase error. If $b$ is a worst-signed additive flag bias, the sufficient active pass condition is

$$
\frac g4-b>0.226073713789.
$$

Useful one-parameter limits are:

| only imperfection | required value |
|---|---:|
| combined scalar contrast | $g>0.904294855157$ |
| two identical lossy sample traversals, otherwise ideal | per-pass power transmission $>0.950944191400$ |
| differential phase only | $|\delta|<0.441071253192$ radians $=25.2715^\circ$ |

The two-pass number is per flag, not compounded over the three sequential photons; the majority formula already accounts for all three flags. Any additional loss, imperfect overlap, phase error, or bias must share the same small budget.

## Feasibility verdict

The $N=4096$ protocol is now resource-complete and has a quantitative active-side tolerance. It remains experimentally ambitious rather than experimentally credible because no platform mapping yet supplies measured $T_L,T_R,\xi,\nu,\delta$, exact $4096$-mode Sylvester transforms, or stable coherent access to the four $4096$-coordinate sign banks.

Nominal addressable-pixel count alone cannot pass this gate. A platform owner must demonstrate the combined contrast inequality with uncertainty and must retain no-clicks in the declared decision record. Retrying until detection would change the hard-dose experiment and is not permitted.

## Reproduction

Run:

```bash
/opt/homebrew/Caskroom/miniconda/base/bin/python3 searches/active_six_robustness_gate.py --write-artifact
```

The committed artifact is `artifacts/active_six_robustness_gate.json`.
