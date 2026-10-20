# Continuous-parameter repair after the finite-bank failure

Chronology: the 21-case common-observable comparison is complete. Its primary
median target RMSE ratio is 1.93, below the frozen twofold gate. The finite-bank
action comparison is also complete: target disagreement costs more than global
disagreement, and all three off-grid target distributions lie outside its
reported envelopes. Those results and gates remain unchanged. This follow-up
repairs a confidence-set limitation; it does not retroactively confirm them.

Replace finite point membership with a continuous supplied domain:
mu in [.03,1.5], pL,pR in [0,1], captured-removal multiplier h in [.15,1.2].
Retain the same single-motor, controlled-source capture family and randomized
snapshot observable. Coverage remains conditional on that family and domain.
Do not claim protection against arbitrary biological misspecification.

Propagate lower/upper acyclic flow intervals over parameter boxes. Use outward
rounding for elementary arithmetic, nonnegativity and probability constraints.
Prune a box only if one of its predicted category intervals is disjoint from
a simultaneous calibration interval. Subdivide retained boxes without dropping
any on resource exhaustion. The union of doubled-removal snapshot enclosures
gives a continuous-domain intervention envelope. This is classical interval
analysis/branch-and-bound applied to the explicit transport model, not a new
general inference architecture. Exhaustive point inclusion checks and nested
box tests must precede experimental evaluation.

New mechanisms, not among the completed 21 cases:
(.21,.3,.7,.4), (.45,.65,.35,.8), (.8,.4,.6,.55); two seeds each.
Acquire 10,000 baseline snapshots at each gamma=.3/.1 and 20,000 independent
target snapshots at gamma=.3 after doubling removal. Per-category calibration
Clopper-Pearson error is .05/(2*82); the union event covers both fixed assays.
Use the same data for the original finite-bank diagnostic, with no off-grid
coverage claim. Compare enclosure/abstention, width, box counts and compute.

Start from one box covering the entire domain. Prune, evaluate target widths,
then split up to 512 of the widest target boxes along their longest normalized
parameter side (ties by parameter/box order). Limit 4096 retained boxes and
32 refinements. On any resource limit return the full remaining envelope and
abstain if its maximum category width exceeds .02. Reject only when no boxes
remain. Stop once every target category width is <=.02. Retain final boxes,
observed counts and all failures. No post-outcome setting changes.

Pass this exploratory repair only if all six generating target distributions
are enclosed and at least four reach .02 precision within the stated budget.
Independently report fresh engine outcomes and Monte Carlo discrepancy. A pass
concerns safe continuous parameter handling, not novel experimental design,
superiority to dLEM, independent measured biology or end-to-end cost savings.
