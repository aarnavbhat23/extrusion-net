# Approximating saturation without a dense contact matrix

This is an experimental specialization of established symmetric tensor lifting
and polynomial kernel approximation, not a new model class. It addresses the
specific gap between the exact linear scan and the saturated contact operator.

## Transfer matrix

For `i < j`, define

```
T_i = [[pR_i, 0], [load_{i+1}, pL_{i+1}]]
D(i,j) = e_1^T T_i ... T_{j-1} e_0.
```

Indices of the two basis vectors are zero-based. Expanding the lower-left entry
of the product recovers the sum over motor load sites in `DERIVATION.md`.

For a positive integer `p`, take the symmetric tensor power of this transfer.
In the monomial basis `q^p, q^(p-1)u, ..., u^p`, its entries are

```
M_p(T)[r,s] = choose(r,s) * a^(p-r) * b^(r-s) * c^s,  0 <= s <= r <= p,
```

where `a=pR_i`, `b=load_{i+1}`, and `c=pL_{i+1}`. This representation respects
matrix products, and

```
D(i,j)^p = e_p^T M_p(T_i) ... M_p(T_{j-1}) e_0.
```

A right-to-left recurrence with injections `e_0 v_j` computes the upper
triangle applied to values. A left-to-right recurrence with transposed transfers
and injections `e_p v_j` computes the lower triangle. The diagonal is excluded.
`contact_power_scan` computes an **entrywise power**, not the ordinary matrix power.
It permits exactly zero pass gates without prefix division.

## Saturation and an explicit limitation

Approximate `1-exp(-x)` with

```
P_q(x) = sum_{p=1..q} (-1)^(p+1) x^p / p!.
```

For `0 <= x <= B`, Taylor's theorem gives the exact-arithmetic uniform error
`B^(q+1)/(q+1)!`. Consequently the error in each message feature is bounded by
that quantity times `sum_j abs(v_j)`; mass error is at most `(n-1)` times it.
These are **absolute**, potentially loose bounds. They do not bound floating-point
roundoff, relative error, or the complete network's task accuracy.

For gates in `[0,1]`, a computable argument bound is

```
B = beta * max(load) * max_{d=1..n-1} d * p_max^(d-1),
p_max = max(max(pL), max(pR)).
```

There are `q(q+3)/2` lifted state coordinates per value feature per direction.
This readable implementation uses dense lifted transfers and costs
`O(n q^3 d_value)` across orders, including increasing constants. It stores
sequence-length arrays and training states; it is not constant-memory training.
For fixed `q`, dependence on sequence length is linear.

Large arguments require high degree and can cause severe cancellation. If all
passes and loads are one, `B=beta*(n-1)`; fixed degree then gives an unusable
bound. The experimental routine intentionally does not replace the trained mixer
or silently clamp contacts. A small-argument numerical success is not permission
to use the approximation for a checkpoint whose argument bound is large.

## Relation to prior work

- [Hydra](https://arxiv.org/abs/2407.09941) and
  [Mamba-2](https://arxiv.org/abs/2405.21060): structured sequence mixers.
- [TaylorShift](https://arxiv.org/abs/2403.02920): polynomial expansion and
  algebraic reordering for efficient attention. Polynomializing an interaction
  kernel is therefore not a first here.
- [hm-toolbox](https://sma.epfl.ch/~anchpcommon/publications/hm-toolbox.pdf):
  hierarchical structured matrices and fast Hadamard products. Structure under
  entrywise products is an established numerical-linear-algebra topic.

The narrow candidate contribution is a specialized, tested recurrence for this
contact parameterization plus an honest account of its useful regime. Priority
and usefulness beyond these established tools are not established by this note.

## Numerical study fixed before its result

Compare orders 1, 2, 4, 6, 8 on seeded nonnegative gates with a bounded argument,
at lengths 32 and 256 in float64. Check exact power messages and gradients against
dense references, observed saturation error against the stated analytic bound,
and timings. Also record an all-open gate case to expose the large-argument
failure. Do not train a new model or select an approximation degree from downstream
test accuracy. No claim of a faster saturated model follows from this study.
