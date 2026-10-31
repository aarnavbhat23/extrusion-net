# Noisy shared-clock recovery and intervention prediction

**Conditional on the geometric shared-clock family, noisy local bands recover
the clock and improve held-out intervention prediction.** At output noise SD
.001, mean clock absolute error is .00198 and mean message RMSE is .00202,
versus .04346 for the static marginal-gate oracle. All 18 mechanisms' clock and
clean held-out message intervals cover in this screen. These are simulated
mechanisms and conservative model-conditional bounds.

[Protocol](CLOCK_IDENTIFICATION_PROTOCOL.md), [inversion](../clock_identification.py),
[all 96 rows](../results/clock_identification_study.json).
The protocol was locally written before outcomes; the screen is exploratory.
Two sizes 16/64, three gate mechanisms each, three shared survival probabilities
and four noise levels supply 72 in-family rows. The nongeometric lifetime mixture
supplies 24 additional mismatch rows. Every empty set and wide interval remains.

## Confidence propagation

Distance-two bands under two known loadings linearly determine effective gates
(a,b)=(s*l,s*r). If every band scalar has error at most eta, each local inverted
gate error is bounded by eta times its row's absolute inverse-design sum.
Intersect gate intervals with [0,1]. For each distance-three observation,
interval arithmetic bounds the static contribution S and middle weight
w=q_middle*a_i*b_(i+1). The equation

    observed_distance_three - S = w * theta,  theta = 1/s - 1 >= 0

gives theta bounds from residual bounds and positive denominator bounds. A zero
lower denominator gives an infinite upper theta bound rather than an unstable
division. Intersect all location/view constraints. Also s>=max(effective-gate
lower bounds), since raw pass gates are at most one. Convert back to s, then
raw gates a/s,b/s. Invalid intersections reject the family.

For N observed scalars with independent Gaussian noise SD sigma, choose
eta=sigma*sqrt(2*log(2*N/.05)); a union bound covers all band errors with at least
95% probability. (Independence is not needed for that union bound if the same
marginal Gaussian bounds hold.) The code adds 1e-14 outward numerical padding.
There are 58 or 250 observed scalars at the two sizes. This is additional
band-resolved access, not the original length-independent vector-query budget.

Positive-contact monotonicity propagates the raw-gate/clock box to clean
held-out messages under new uniform loading, a left cut and simultaneous left/
right cuts. The same calibration event covers every such prescribed prediction.
Raw-gate/clock correlations are discarded, so the resulting prediction intervals
can be much wider than necessary. No noisy-test-output coverage is claimed.

| Band noise SD | Clock MAE | Clock-model message RMSE | Static oracle RMSE | Mean prediction width | In-family covered / total | Mixture rejected / total |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | below 2e-15 | below 4e-16 | .04346 | below 1e-12 | 18/18 | 6/6 |
| .0001 | .000198 | .000202 | .04346 | .03275 | 18/18 | 6/6 |
| .001 | .001977 | .002018 | .04346 | .33095 | 18/18 | 6/6 |
| .01 | .019933 | .020352 | .04346 | 10.16069 | 18/18 | 0/6 |

The static oracle knows the exact noiseless effective marginal gates. It is a
stronger comparator than fitting them from noisy bands, and its informational
advantage is explicit. The point estimate uses a median ratio; its empirical
error is separate from the conservative confidence box. Coverage here is a
sanity check over 18 mechanisms per noise level, not a powered calibration study.
Bounded-noise numerical tests also check clock/raw-gate enclosure.

## Decision and limits

The point-recovery payoff gate passes. Continue the bounded identification
problem, while stopping the unsuccessful causal clock learning variant described
in [its separate report](CAUSAL_CLOCK_RESULTS.md). At the largest tested noise,
prediction bounds become impractically wide and miss every mixture mismatch;
that limitation remains central. Nongeometric lifetimes, heterogeneous speed,
unknown loading and measured chromatin still require a different model set or
new evidence. The next useful theory is a tighter joint feasible set for
intervention predictions, compared at identical band-observation cost, followed
by transfer to an independently implemented dynamic engine. This result does
not establish a generally superior sequence learner or new biological physics.
