# Continuous parameters: valid enclosure, failed precision and feasible ambiguity

**Continuous-domain enclosures retain all six generating mechanisms, but none
meets the .02 maximum-category-width gate.** A survival identity tightens the
bounds without solving precision. Numerical witnesses, independently checked
with 70-digit arithmetic, show why computation alone is insufficient: in every
case, two mechanisms satisfy all chosen calibration inequalities yet differ by
.0302–.0474 in intervention absence probability.

This is an obstruction for envelopes covering every parameter consistent with
the chosen rectangular calibration constraints. It is not an information-theoretic
impossibility for every statistical procedure, nor a biological discovery.

## Chronology and retained failures

The [continuous-domain protocol](CONTINUOUS_OCCUPANCY_PROTOCOL.md) was committed
at `833a651` after the [finite-bank failure](OCCUPANCY_RESULTS.md), before these
new cases. The interval implementation at `0e93778` preceded evaluation.
All six initial runs enclosed the generating parameters and intervention
probabilities but abstained, failing the requirement that at least four achieve
.02 width. They remain in [the original results](../results/continuous_occupancy_study.json),
with final parameter boxes in [the array archive](../results/continuous_occupancy_boxes.npz).

The [aggregate-survival amendment](CONTINUOUS_AGGREGATE_AMENDMENT.md), committed
at `d1cf3e6` after seeing that absence dominated every envelope, specifies a
theory-driven reanalysis. This is the same six datasets, not fresh confirmation.
[Refined results](../results/continuous_occupancy_refinement.json) and
[their parameter boxes](../results/continuous_occupancy_refined_boxes.npz) are
saved separately. The frozen settings and precision target are unchanged.

## Continuous-domain enclosure

The supplied domain is mu in [.03,1.5], pL,pR in [0,1] and h in [.15,1.2].
The assumed process is still one motor, known loading, convergent barriers,
unit step rates, and captured-state-dependent removal. Arbitrary biological
misspecification, traffic and measurement biases are not covered.

Each new mechanism uses two independent 10,000-episode calibration assays at
gamma=.3/.1, and a fresh 20,000-episode target at gamma=.3 after doubling removal.
Clopper-Pearson intervals use error .05/(2*82) for each categorical marginal.
Dependence among categories does not invalidate the union bound. Conditional
on the family/domain and independent trials, simultaneous calibration enclosure
is at least .95 per case, not across the full six-case study.

Outward-rounded elementary arithmetic propagates nonnegative flow intervals
over each parameter rectangle. A box is pruned only when some predicted
category interval is disjoint from its observed calibration interval. Splitting
retained boxes preserves their union; hitting the 4096-box/32-refinement cap
never licenses dropping a remaining box. The intervention envelope covers every
retained box. Thus a generating parameter remains on the simultaneous calibration
event, and its target remains enclosed. These are interval-analysis and
set-membership methods, not a new generic inference architecture.

## Aggregate survival identity and exact rectangle corners

Let S be the first successful capture time, infinite if neither barrier captures.
Before S removal has rate mu; afterward it has rate mu*h. With s=mu+gamma,

    K = E[exp(-s*S); S finite]
    P(present) = gamma/s*(1-K) + gamma/(mu*h+gamma)*K.

The first term is observation before capture; the second is capture before
observation/removal followed by survival to the residual observation clock.
Memorylessness justifies the residual competing-clock calculation.

For independent unit-rate arms and barriers three steps away,

    Ksingle = 1/(s+1)^3
    Kboth = 2/(s+2)^3 + 6/(s+2)^4 + 12/(s+2)^5
    K = [pL*(1-pR)+pR*(1-pL)]*Ksingle + pL*pR*Kboth.

Kboth is the Laplace transform of the minimum of two Erlang-three clocks.
Equivalently, K=pL*H(3,0;s)+pR*H(0,3;s)-pL*pR*H(3,3;s), connecting the earlier
reach calculation to snapshot occupancy by inclusion-exclusion. These are
classical competing-risk/Laplace-transform identities, not new extrusion physics.

Couple two processes to the same step/capture/observation randomness. Increasing
mu or h increases cumulative removal hazard, so survival decreases. Increasing
a capture probability makes first capture earlier. This increases survival for
h<=1 and decreases it for h>=1. Therefore rectangular extrema of total absence
occur at the prescribed mu/h endpoints and the corresponding capture corners,
including boxes crossing h=1. Outward evaluation at those corners tightly bounds
the aggregate without summing 81 dependent state intervals. Tests compare this
identity with the independent full-state flow and check sampled box inclusion
on both sides of h=1.

## All six outcomes

Widths and attained spans below are percentage points. The precision target is
2 points in every category. The last column counts boxes retained after the
aggregate refinement, not point estimates.

| mu / pL / pR / h | Seed repeat | Initial max width | Refined max width | Feasible absence span | Refined boxes |
|---|---:|---:|---:|---:|---:|
| .21 / .3 / .7 / .4 | 0 | 18.72 | 14.50 | 4.50 | 2008 |
| .21 / .3 / .7 / .4 | 1 | 21.50 | 17.77 | 4.74 | 1742 |
| .45 / .65 / .35 / .8 | 0 | 8.48 | 7.55 | 3.99 | 3839 |
| .45 / .65 / .35 / .8 | 1 | 10.11 | 7.62 | 4.19 | 4096 |
| .8 / .4 / .6 / .55 | 0 | 5.59 | 5.21 | 3.13 | 4096 |
| .8 / .4 / .6 / .55 | 1 | 5.63 | 5.15 | 3.02 | 4096 |

The refinement improves maximum width by a median factor of 1.17 and retains
all generating parameters/targets, but 0/6 achieve useful precision. Remaining
widths exceed feasible witness spans, so computational conservatism also remains.
All original finite-bank candidates reject on these two-assay datasets; that
is different from the earlier bank policy's premature precision stops.

## Post-outcome ambiguity diagnostic

After the refinement, a constrained numerical search minimizes/maximizes the
intervention's absent probability within the same domain and calibration
inequalities. It uses three starts, including the generating parameter solely
as a diagnostic feasible witness. This knowledge is not supplied to an estimator
or acquisition policy. Saved pairs have nonnegative calibration slack; no claim
of globally optimal extrema follows from the optimizer.

An independently implemented Decimal state-flow calculation with 70-digit
precision verifies all 12 supplied witness populations against every represented
calibration inequality and verifies all six intervention spans exceed .02.
[Witness populations](../results/occupancy_ambiguity_diagnostic.json) and
[independent arithmetic checks](../results/occupancy_witness_verification.json)
are retained. Any envelope containing every parameter feasible under these
rectangular constraints must contain both witnesses. More subdivision cannot
meet .02 on these datasets without changing the constraints or observation model.

This does not show that the datasets are insufficient for every valid statistical
method. A joint likelihood-based confidence construction might exclude points
that conservative marginal rectangles retain. Direct capture/residence-time
measurements could also resolve a different ambiguity. Test these possibilities
before adding model capacity or scaling a predictor.

```bash
python study_continuous_occupancy.py
python refine_continuous_occupancy.py
python diagnose_occupancy_ambiguity.py
python verify_occupancy_witnesses.py
python plot_occupancy_research.py
```
