# Paper decision after the finite-data calibration round

**2026-10-29 continuation:** full manuscript work remains deferred. The
[competing-loss direction](LOSS_PERTURBATION_RESULTS.md) passes exact-data
feasibility, but does not yet change the grades below. Keep theory notes as
research records; prioritize finite-data and observation-model validation.

Date: 2026-10-28. **Write a short internal theory draft; do not submit yet.**
Do not frame this as a replacement for neural networks or attention. Current
results support a precise question about information lost by an observation
protocol, with a useful negative finite-data result.

Working title: **What snapshot normalization destroys in loop-extrusion inference**.

The central claim is that a known exponential observation clock does not fix
the scale of a killed kinetic process when only normalized retained-state
histograms are observed and motor speeds are unknown. An explicit family of
observationally equivalent models predicts different removal interventions.
For the specified capture model, two clocks or one absolute-retention snapshot
restore exact identification. At the tested finite counts, the simpler
absolute-retention observation contract gives substantially smaller point errors.

This is an observation-contract result. The elementary resolvent transformation,
two-point statistical argument and binomial confidence construction are not
claims of new general mathematical tools. Priority for their exact combination
and application remains unestablished.

## Draft outline and claim ledger

1. Define known loading, random-time sampling, hidden capture and the information
   actually retained. Contrast this with recorded initiation counts and timestamps.
2. State and prove normalized-resolvent ambiguity. Include restrictions on allowed
   parameter domains: the transformed mechanism must remain inside the domain.
3. Give intervention separation and standard two-point lower bounds.
4. Derive the capture-family inverse for two clocks and one absolute denominator;
   handle zero capture and irrelevant h explicitly.
5. Report all 72 fits, the failed primary gate, unbounded confidence sets and
   actual information/cost differences. Include the strong full-retention control.
6. State model/observable limitations and what a real experiment would require.

The existing conditional source-recovery theorem is a separate result. It may
be referenced for context but should not be combined with many unrelated failed
architectures into one unfocused manuscript. The strongest paper needs one
question, one observation contract and one defensible conclusion.

## Requirements before submission

- Complete equation-level priority review against kinetic identifiability,
  killed-process conditioning and snapshot-calibration results. Obtain external
  mathematical feedback. A search that fails to locate the formula proves nothing.
- Compare timestamp likelihood, absolute-retention calibration and standard
  information-based design using the same observable and an explicit cost ledger.
  If retention still wins, retain that answer instead of inventing a new variant.
- Address uncertain loading and state-dependent dropout. Derive a restricted
  theorem or demonstrate failure under each; never imply the exact inverse survives
  untested violations. Freeze new regimes before results.
- For a useful finite-data method claim, require a prespecified intervention
  error tolerance and coverage check on fresh mechanisms, not just optimizer
  convergence. The present speed intervals do not meet that requirement.
- For a biological paper, validate the distance-to-anchor observable and use a
  genuinely withheld perturbation. Existing simulator transfer does not suffice.
  A theory note can proceed without that claim if its theorem is substantial and
  original enough after review.

## Current assessment

Subjective evidence-based assessment, not a peer-review prediction: novelty
about **6.5/10**, overall research quality **7.5/10**, paper readiness **6.5/10**.
Proof/test discipline is around 8 and reproducibility around 9; independent
biological validation is only around 4–5. These are not an arithmetic average:
novelty and applicability are bottlenecks. Uncertainty in these judgments is
at least roughly one point. Failed gates improve credibility but do not create
a new capability or justify a higher novelty rating.

Relevant primary precedents include [kinetic identifiability and experiment
selection](https://academic.oup.com/bioinformatics/article/35/24/5216/5514038)
and [time-scale ambiguity in stationary spatial snapshots](https://arxiv.org/abs/2610.07341).
The latter studies reaction-diffusion spectra, a different observation model.
The comparison narrows potential originality; it does not establish it.
