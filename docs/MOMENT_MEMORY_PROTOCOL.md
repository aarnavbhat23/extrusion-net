# Moment-conditioned retrieval: equations, novelty boundary and test gates

Research hypothesis: can a fixed-size positive-feature memory retrieve an atom
and expose ambiguous answers more usefully than weighted averaging? This is an
operator-level candidate, not yet an end-to-end sequence architecture. Key/value
bindings are supplied to the operator. Learned binding from a raw token stream
is an additional gate. No claim of priority or a 9/10 result is made.

## Equations and existing ingredients

Let nonnegative unit-normalized key features be phi_theta(k), and let
xi(v)=(1,v). Store the symmetric tensor

    S_t = S_(t-1) + phi_theta(k_t) tensor xi(v_t) xi(v_t)^T.

At query q, contract G(q)=sum_a phi_theta(q)_a S_a + lambda I. For a catalogue
candidate value v, read its score d(v)=xi(v)^T G(q)^(-1) xi(v). Point decoding
chooses the smallest score. If q exactly equals a written key and its write
weight is one, that associated value necessarily satisfies G-xi(v)xi(v)^T>=0,
equivalently d(v)<=1 for positive-definite G. Return a certified value only
when exactly one catalogue value satisfies this test; otherwise abstain.

This is a conditional certificate. Absent keys, inexact query/key alignment,
unknown write weights, and predicting a future label are not covered. Duplicate
keys with conflicting values must remain ambiguous. Ridge regularization adds
PSD mass and can only loosen the necessary condition. Numerical tests need
outward tolerance; exact-real arithmetic supplies the mathematical statement.

The tensor update is an additive lifted-feature memory. The inverse-moment score
is an empirical Christoffel/Mahalanobis construction. Neither is a new general
mathematical primitive. The candidate distinction is the query-conditioned value
moment, atom-compatibility read and its conditional retrieval use. A targeted
search found no exact matching formulation; this is not an exhaustive priority
search. Further literature/code comparison remains mandatory before a novelty claim.

Closest precedents:
- [HOLA](https://arxiv.org/abs/2607.02303): residual-selected exact cache. This
  rules out claiming a residual cache itself as novel.
- [Sparse Delta Memory](https://arxiv.org/abs/2607.07386): sparse expansion of
  recurrent memory. Expanded capacity alone is not the proposed distinction.
- [Memory Caching](https://arxiv.org/abs/2602.24281): growing recurrent memory.
- [Test-time regression](https://arxiv.org/abs/2501.12352) and
  [Variational Linear Attention](https://arxiv.org/abs/2605.11196): regression
  and inverse-state approaches. Rebranding RLS is excluded.
- [Empirical Christoffel functions](https://arxiv.org/abs/1701.02886): direct
  precedent for inverse-moment support scores. The atom PSD identity is elementary.
- [Higher-order linear attention](https://arxiv.org/abs/2510.27258) and
  [curvature-conditioned queries](https://arxiv.org/abs/2606.01294): second-order
  statistics already have sequence-model precedents; inspect distinctions.

## Frozen operator screen

Raw keys have 16 continuous coordinates. Train, validation and test key pools
are generated independently and are disjoint. Values are the eight 3-bit codes.
Train on eight bindings per episode, batch 16, 600 Adam steps, seeds 0/1/2.
Validation is in-distribution in number of bindings but uses unseen raw keys.
Evaluate 8/16/32/64 bindings, 64 episodes per size, without test-based selection.

All learned encoders have a common 16-by-16 transform and a temperature scalar;
fixed random positive address expansions vary to match actual stored-state
budgets. Candidate address width 32 stores 32*4*4=512 scalars (full matrices).
Controls: normalized averaging (128*4=512), delta-rule memory (170*3=510),
RLS (21*21+21*3=504), full softmax attention and a FIFO cache capped at the same
512-scalar budget (26 raw 16-coordinate keys plus 3-coordinate values). Count
projection buffers, learned parameters and stored state separately. Also retain
same-address-width averaging as a mechanism ablation, with its smaller state.

Record point accuracy, certified coverage, accuracy counting abstention as failure,
false certificates, update/read time, and state. Include absent queries and
conflicting duplicate writes as assumption violations. Do not label absent-key
certifications as valid, or interpret ambiguity as successful full retrieval.
Continue the architecture claim only if mean length-32 accuracy is >=95%,
certified coverage >=90%, no false certificates occur under the valid assumptions,
and matched-state controls do not dominate quality/resources. A pass needs fresh
confirmation and learned token-stream binding, not just this supplied-binding test.

## Small real-text diagnostic, distinct from confirmation

A public Shakespeare character corpus can establish a causal text harness and
measure practical costs. Freeze a contiguous 80/10/10 split before training,
training-only vocabulary, block length 64, batch 8, 600 steps and seeds 0/1/2.
Compare the moment read, matched-state averaging, static transport and causal
attention on a shared small backbone. Select only by validation next-character
cross entropy, report untouched test cross entropy, parameter counts and timing.
Contextual moment writes may use only already observed tokens; test prefix
invariance. This tiny diagnostic is not a language-model competitiveness result
and cannot rescue a failed operator gate. No atom certificate applies to the
correctness of a future character. No scaling follows from this pilot alone.
