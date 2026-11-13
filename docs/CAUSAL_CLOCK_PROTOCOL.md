# Causal shared-clock model screen

Exploratory architectural screen. Written before training outcomes. Keep the
original released models and results intact. The candidate adds a capped shared
geometric lifetime correction to a global causal static transport. A matched
independent-arm correction has the same lifetime parameter, age cap and learned
correction strength. A static transport and causal softmax attention are controls.
All use the same one-block causal convolution/MLP backbone, width 32 and four heads.
The age cap is 8, restricting the correction to at most 17 positions; the static
path remains global. This is a finite-range model extension, not an implementation
of unbounded-rank infinite-lifetime transport or a general memory solution.

Train on the existing recall and bracket generators at lengths 24/32/48/64,
800 sequences, batch 24, AdamW lr .002 and 400 steps. Use seeds 0/1, the same
batches per seed across models and in-distribution validation every 100 steps.
Select checkpoints only by validation. Evaluate lengths 64/128/256, 128 sequences
each. Retain all model/task/seed results and timing. No hyperparameter tuning from
long-length outcomes. This pilot is too small to establish language-model value.

Continue this architecture only if it beats both static transport and the matched
independent-clock control by at least two percentage points in mean length-128/256
accuracy on at least one task, without more than two points in-distribution loss.
A positive result needs a separately frozen confirmation with more seeds, stronger
controls and real data. A negative result should stop this version, not trigger
unreported tuning. Count parameters, throughput and the explicit age-dependent cost.
