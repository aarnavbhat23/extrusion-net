# A causal clock correction does not solve the model's learning limitation

**Decision: stop this architectural variant as an improvement claim.** The
shared-clock correction failed its prewritten comparison rule against both a
static causal transport and a matched independent-clock correction. All 16
training runs, including the failures, are retained.

[Protocol](CAUSAL_CLOCK_PROTOCOL.md), [implementation](../causal_transport.py),
[training results](../results/causal_clock_study.json),
[resource screen](../results/causal_transport_resources.json).
This is a small two-seed synthetic screen. The protocol was written locally
before training, not independently committed as a prospective confirmation.

## Actual model change

The new backbone is causal throughout: left-padded local convolution, causal
transport and no future-looking operation. The shared clock adds a mixture of
capped-age two-stage transport kernels to a global static transport. Its survival
and correction strength are learned per head. A matched independent-clock
correction has the same parameter count and reach cap. The age cap 8 limits the
correction to 17 positions; the static path remains global. It is not an
implementation of the infinite shared-clock family's arbitrary-rank kernel.

An associative-doubling implementation evaluates gated prefix/window sums
without building an n-by-n transport matrix. For fixed width and head count,
the current static implementation takes O(n log n) work; the capped clock adds
O(n H log H) work for age cap H. This is not an optimized linear-time kernel.
Autograd retains intermediate arrays. Parameters, activations and finite-age
costs must be counted separately from the original sequential scan.

Tests compare kinetic messages against the independent contact enumerator,
check clock gradients by finite differences, and verify prefix invariance for
the complete model when all future tokens change.

## Results

Mean accuracy over two seeds; models were selected using only short-length
validation. Training lengths were 24/32/48/64, with 400 optimizer steps and
800 training sequences. Attention has about 1,580 more parameters in this small
backbone, so this is a reference control, not an exactly equal-parameter trial.
The three transport variants have the same parameter count; the static variant's
clock/correction parameters are unused. This tiny difference cannot explain the
near-total retrieval gap but should not be hidden.

| Task/model | Length 64 | Length 128 | Length 256 |
|---|---:|---:|---:|
| Recall, static | .3083 | .2915 | .2316 |
| Recall, shared clock | .3098 | .2930 | .2332 |
| Recall, independent clock | .3129 | .2930 | .2332 |
| Recall, causal attention | .9969 | 1.0000 | .9968 |
| Brackets, static | .9553 | .9193 | .9051 |
| Brackets, shared clock | .9519 | .9101 | .8925 |
| Brackets, independent clock | .9499 | .9136 | .8969 |
| Brackets, causal attention | .9427 | .9017 | .8585 |

Static transport does better on bracket extrapolation in this screen. Shared
clock is slower and supplies no useful retrieval improvement. This does not
show that all lifetime models are useless; it rejects the particular candidate
and training regime tested. No long-length tuning was performed after seeing
these outcomes.

The post-screen resource check times forward/backward at lengths 64/256/1024/4096
on CPU with two threads and batch two. It records median tokens/s and the total
number of saved nonparameter tensor entries. Those entries can double-count
aliases and are not peak allocated memory. Resource results are implementation-
and hardware-specific and were not used for model/checkpoint selection.

At length 4096, measured training throughput is approximately 331k/119k/230k/93k
tokens/s for static/shared/independent/attention respectively. Shared clock
retains about 38.8 million nonparameter tensor entries versus 4.9 million for
the fused attention implementation. Faster throughput at that length does not
compensate for its accuracy failure or establish a memory-efficiency advantage.

## Direction toward a stronger architecture

The current limiting capability is content-dependent memory selection. More
lifetime parameters or an enlarged finite context do not solve it in these
experiments. A broadly useful model must demonstrate content retrieval,
long-context retention, causal execution and a favorable quality/resource
tradeoff on independent tasks and real data. An impressive name or a biological
analogy establishes none of those properties.

Selective state updates and associative memory already have major precedents:
[Mamba](https://arxiv.org/abs/2312.00752),
[DeltaNet](https://arxiv.org/abs/2406.06484) and
[Titans](https://proceedings.neurips.cc/paper_files/paper/2025/file/a4ca07aa108036f80cbb5b82285fd4b1-Paper-Conference.pdf).
Adding a content gate, a delta-rule write or an attention path would therefore
be a necessary strong control, not automatically a new contribution.

The next defensible architecture hypothesis is whether **memory capacity can
expand only when a transport state cannot represent conflicting associations**,
with a testable cost/retention rule and controlled-probe audit. This is a research
question, not a specified new model or a novelty claim. Before implementation,
write the update rule, its relation to dynamic memory/data structures and
existing memory models, and an exact small-task oracle. Compare against simple
key-value storage, selective state models and content attention under the same
state and compute budget. Reject the direction if it reduces to an existing
method or wins only when supplied the test keys. Any positive tiny-task result
would still require language-model and real-data evaluation before scaling.
