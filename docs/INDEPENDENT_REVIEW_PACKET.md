# Review packet: predicting a controlled-loss target under loading uncertainty

Status, 2026-11-01: prepared for outside review; **no external mathematical,
priority or experimental feedback has been obtained or requested by this
repository**. Computational replay is not external peer review. Do not submit
biological-identification claims on the current evidence.

2026-11-03 follow-up: an [exact-arithmetic audit and measurement contract](MEASUREMENT_CONTRACT_AUDIT.md)
confirms the counterexample and the sharpness proof within scope, corrects three scope statements,
shows that the current calibration contract cannot beat naive dose extrapolation at any count, and
freezes go/no-go criteria for a reconstitution feasibility pilot. No pilot data exist.

## Proposed paper core

Question: when do retained/initiated observations and a separately calibrated
reference bound one intervention response even though individual kinetic rates
and loading are unidentified?

Read only these core artifacts first:

1. [Loading-target theory](LOADING_TARGET_THEORY.md): the complete-snapshot gauge,
   scalar target identity, sharp loading-drift polytope and endpoint mechanisms.
2. [Counted-reference theory](COUNTED_LOADING_THEORY.md): four-probability
   confidence image, sensitivity counterexample, classical-equation equivalence.
3. [Original finite-noise results](LOADING_RESEARCH_RESULTS.md) and
   [counted-reference results](COUNTED_LOADING_RESULTS.md), including failed
   same-loading assumptions and timestamp refinement.
4. `verify_loading_research.py` and `verify_counted_loading.py`: optimizer-free
   checks of the recorded primal/dual certificates and endpoint mechanisms.

The original attention replacement and passive-memory explorations remain
separate. These results do not establish a new learning architecture.

## Questions for a mathematical reviewer

* Does the fixed-gamma gauge preserve the entire claimed emission/absence law,
  including the role of loading and optical loss? Which variants violate it?
* Are the reduced-polytope necessity AND endpoint sufficiency proofs complete,
  including W=0, nonnegative rates, ell=ell_max and kappa<=1? Do not silently
  extend them to kappa>1, unknown target loading or h!=1.
* Is the confidence-image guarantee stated at its true level (a conservative
  rectangle image), rather than an optimal joint confidence set?
* Is the counted-reference counterexample a valid complete observation-law
  equivalence under its Bernoulli model?

## Questions for a priority reviewer

The reciprocal response curve is algebraically Stern–Volmer shaped; it is not
the originality claim. Partial identification through LPs is also established.
Does existing fluorescence/kinetic experimental-design work already give this
specific controlled-loss target, loading-drift polytope or equivalent protocol?
Which earlier result is the closest? What remains worth publishing after that
comparison? A literature search that finds no exact title match is insufficient.

## Biological eligibility and a feasible pilot decision

A cohort is eligible only if initiated denominators, absolute retention in both
arms, an independently checked intervention factor, independent time-zero
loading-reference counts, and reference sensitivity/false-positive information
are recorded. Cells/movies, not individual correlated spots, must define the
independent units unless a justified hierarchical observation model is used.
The target cohort must preserve baseline loading and the observation operator.

Existing normalized iFRAP curves and locus-distance trajectories do not supply
these measurements. GRID's public CDX2 metadata concerns illumination-dependent
survival; it does not supply this intervention/loading contract. Its archive was
not downloaded (website 403; official file API 401), so it contributes no new
calibration outcomes here. No source data were fabricated or substituted.

Before collecting a large cohort, an imaging collaborator should judge whether
an initiated loading episode and pre-loss reference can actually be observed
for this biological system. A controlled reconstitution may make denominators
more concrete, but applicability, independent reporter specificity, the
state-independent h=1 hazard and a calibrated linear perturbation must first be
demonstrated. The current event-count grid is a statistical design calculation,
not a claim that 64,000 independent reference events are experimentally cheap.

Freeze specimen allocation and treatment labels before target-outcome analysis;
fit reference/error operators on separate specimens; blind the intervention
response until prediction intervals are locked. If eligibility fails, keep the
paper theoretical and explicitly illustrative. Do not impute the missing
loading measurements from the same survival observations they must calibrate.

## Submission decision

The packet supports a focused theory/methods discussion, with substantial
computational checks and retained negative findings. External priority/proof
feedback remains open. Biological-response identification is not validated.
Choose a venue and manuscript scope only after the contribution survives that
feedback; do not present project grades as acceptance probabilities.
