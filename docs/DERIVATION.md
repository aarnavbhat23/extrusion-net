# Derivation

## Process

Positions are `0..n-1`. A motor loads on the bond just left of position `k` with propensity `l_k`. Its left arm then visits `k-1, k-2, ...`; at each position `m` it stalls with probability `sL_m`. Its right arm visits `k, k+1, ...` and stalls with probability `sR_m`. The arms move independently.

Positions `i < j` are held in the same loop when the load bond lies between them and each arm reaches its end of the pair:

```
D(i, j) = sum_{k=i+1..j} l_k * prod_{m=i+1..k-1} (1 - sL_m) * prod_{m=k..j-1} (1 - sR_m)
```

The arm reaching `i` does not need to pass `i` itself; stalling at `i` still holds `i`. `D` is symmetric by definition and the diagonal is unused.

Many motors load over time. Treating arrivals as Poisson with rate proportional to `D`, the chance that at least one motor has held the pair is `W(i, j) = 1 - exp(-beta * D(i, j))`, with a learned `beta > 0` per motor type. Position `i` receives `sum_j W(i, j) v_j` and `log(1 + sum_j W(i, j))`.

## O(n^2) evaluation

Let `a_m = log(1 - sL_m)`, `b_m = log(1 - sR_m)`, `A_t = sum_{m<=t} a_m`, `B_t = sum_{m<=t} b_m`, and `A_{-1} = B_{-1} = 0`. Then

```
prod_{m=i+1..k-1} (1 - sL_m) = exp(A_{k-1} - A_i)        (k >= i+1)
prod_{m=k..j-1}   (1 - sR_m) = exp(B_{j-1} - B_{k-1})    (k <= j)
```

so

```
log D(i, j) = B_{j-1} - A_i + logsumexp_{k=i+1..j} c_k,    c_k = log l_k + A_{k-1} - B_{k-1}
```

For every row `i`, the range `logsumexp` over `k = i+1..j` for all `j` is one `logcumsumexp` over the masked row. That gives every pair in O(n^2) time and memory, the same order as attention. `tests/test_models.py` checks this against the direct triple sum.

Numerics: stall logits are clamped at 12 (a barrier passes at most about 6e-6 of arms), so `|A_t|, |B_t| <= 12n`. Masked entries use `-1e9` rather than `-inf`, because `-inf` produces NaN gradients in `logcumsumexp`. At `n = 512`, float32 cancellation in `B_{j-1} - A_i + R` is about 4e-4 in log space, a 0.04% relative error in `D`.

## Exactness for the convergent-CTCF rule

**Claim.** Let `sL_m = 1` exactly at forward sites (`>`) and 0 elsewhere, `sR_m = 1` exactly at reverse sites (`<`) and 0 elsewhere, and `l_k > 0` everywhere. Then `D(i, j) > 0` if and only if no `>` precedes a `<` strictly between `i` and `j`.

**Proof.** A term `k` is nonzero iff there is no `>` in `[i+1, k-1]` and no `<` in `[k, j-1]`. Let `f` be the first `>` strictly between `i` and `j` (or `j` if none) and `r` the last `<` strictly between them (or `i` if none). The term is nonzero iff `k <= f` and `k >= r + 1`, so some valid `k` in `[i+1, j]` exists iff `r < f`. That is, every `<` between them comes before every `>` between them, which is the same as saying no `>` precedes a `<`. ∎

So with hard barriers and uniform loading, one extrusion layer computes exactly the insulation pattern of the convergent rule. The magnitude of `D` still depends on spacing, and `W = 1 - exp(-beta D)` approaches an indicator as `beta` grows. `tests/test_models.py` checks the claim on 200 random sequences.

## Structured-matrix view and an O(n) algorithm

Write `L[i, k] = prod_{m=i+1..k-1} (1 - sL_m)` for `k > i` and `R[k, j] = prod_{m=k..j-1} (1 - sR_m)` for `j >= k`. Both are one-sided transports with cumulative-product decay, which makes them 1-semiseparable matrices. The upper triangle of the contact matrix is then

```
D_upper = L diag(l) R
```

the product of two one-sided transports through the diagonal load gate. Any block of `D` whose rows all lie left of its columns has rank at most 2: split the load sum at the block boundary, and each part factors into a function of `i` times a function of `j`. So `D` is 2-quasiseparable. `tests/test_models.py` checks the rank numerically.

This places the *linear* part of extrusion in a known family. Quasiseparable and semiseparable matrices are the structured matrices behind bidirectional state-space models (Hydra, Hwang et al., NeurIPS 2024) and Mamba-2's state-space duality (Dao and Gu, ICML 2024). The one-arm ablation is close to a single bidirectional selective scan. What is specific to extrusion is the parameterization:
- the same stall sequences serve both arms, crossed between the two directions;
- a load gate sits between the two scans;
- the matrix is symmetric by construction;
- the Poisson saturation `1 - exp(-beta D)` is applied elementwise, which leaves the structured family.

It also gives a linear-time algorithm for the unsaturated contact. With `pR_m = 1 - sR_m` and `pL_m = 1 - sL_m`:

```
upper part, right to left:
  Q_k = v_k + pR_k Q_{k+1}
  U_i = l_{i+1} Q_{i+1} + pL_{i+1} U_{i+1}

lower part, left to right:
  G_k = v_{k-1} + pL_{k-1} G_{k-1}
  H_i = l_i G_i + pR_{i-1} H_{i-1}

sum_j D(i, j) v_j = U_i + H_i
```

Four multiply-add recurrences with factors in `[0, 1]`, so they are numerically stable without log space. `extrusion_linear_scan` implements them, and the tests check it against the dense matrix. `bench_scaling.py` shows the scan handles 65,536 positions on one CPU thread, where the dense matrix would need 16 GB (`results/scaling.json`). The `extrusion_linear` mixer trains with `W = beta D` (no saturation), so it can be computed this way. Whether dropping saturation costs accuracy is an experiment in `docs/PROTOCOL.md`.

## The one-arm ablation

Arms start at the query position itself, with no separate load site:

```
W1(i, j) = prod_{m=i+1..j-1} (1 - sR_m)   for j > i
W1(i, j) = prod_{m=j+1..i-1} (1 - sL_m)   for j < i
```

A single such layer cannot represent the convergent rule. Whether a `<` blocks the path from `i` to `j` depends on whether a `>` lies between `i` and that `<`, but a per-position stall probability can't depend on where the query is. Depth can work around this: later layers see features computed by earlier ones. That is what the 1-layer vs 3-layer results show. The form is close to stick-breaking attention with per-position rather than per-pair break probabilities.
