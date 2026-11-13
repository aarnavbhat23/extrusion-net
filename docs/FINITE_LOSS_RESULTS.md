# Finite-data prediction with competing losses

2026-10-29. **The normalized three-condition design fails its prediction gate.**
Retaining absences and scheduled time-resolved snapshots both pass in the matched
model. Exact identification did not imply useful precision for the proposed
normalized design. Full manuscript work remains deferred.

The [protocol](FINITE_LOSS_PROTOCOL.md) was committed as `e4ebb15` before fitting.
All 90 fits and 270 optimizer starts are retained in
[`finite_loss.json`](../results/finite_loss.json). Three new seven-parameter
mechanisms, two episode budgets and fixed seeds were used. Training and holdout
counts are independent; generating parameters are never supplied to fitting.

## Results

Median absolute error in intervention absence, in **percentage points**:

| Setting | Normalized three conditions | Retained absences | Scheduled timestamps |
|---|---:|---:|---:|
| Matched, 6,000 initiated episodes | 5.184 | 1.001 | 0.556 |
| Matched, 24,000 initiated episodes | 2.501 | 0.312 | 0.298 |
| Actual calibration factor 1.8, assumed 2; 24,000 | 2.755 | 1.495 | 1.590 |
| State-dependent extra loss; 24,000 | 1.429 | 0.410 | 0.132 |

At 24,000 matched episodes the per-mechanism medians meet the two-point
tolerance in **1/3, 3/3 and 3/3** mechanisms, respectively. The frozen gate
requires at least 2/3. Under inaccurate perturbation strength these counts are
1/3, 2/3 and 2/3. The stresses have only two repeats per mechanism and cannot
establish broad robustness. Their small errors in some cases are not evidence
that the model mismatch is harmless for other interventions.

All optimizer starts reported convergence. This does not prove a global optimum.
The independent full-state replay checks all counts, all start likelihoods,
holdout scores, predictions and summaries. Maximum likelihood discrepancy is
5.46e-12. The fitting code uses a 100-state reachable generator with nonnegative
uniformization, while data and replay use a separately assembled 324-state
generator with a resolvent or matrix exponential. Agreement verifies computation
within the same mathematical model, not independent biological validity.

## Information and cost

The three methods have different observation contracts. At 24,000 episodes,
half are training and half holdout. Expected summed observation durations are
133,333, 80,000 and 45,000 model-time units, respectively. Scheduled observation
uses eight conditions, versus three and two. This is neither measured scanner
time nor a complete experimental cost model. Exact scheduled timestamps are
not individual random exponential timestamps. No cost-optimal design is proved.

At a *fixed* condition, let retention probability be r(theta) and the
conditional retained distribution be q_i(theta). The absolute distribution is
p=(r q_1,...,r q_d,1-r). For N independent initiated observations with M retained,
the likelihood factorizes into a binomial retention term and a conditional
multinomial term. Under standard differentiability assumptions, per-episode
Fisher information decomposes as

    I_absolute = r I_conditional + grad(r) grad(r)^T / [r(1-r)].

The cross term vanishes because the conditional score has expectation zero.
Thus dropping retention removes a positive-semidefinite information term at the
same sampling schedule. The compared designs have different schedules, so this
identity alone does not rank them. This is an elementary information
decomposition, not a novel theorem. It explains why another normalized design
needs a practical reason to discard the denominator.

## What the confidence check does and does not show

The holdout set is C={theta: L_hold(theta) >= .05 L_hold(theta_hat_train)}.
Conditional on training, the likelihood ratio L_hold(theta_hat_train)/L_hold(theta0)
has expectation at most one under a correctly specified theta0. Markov's
inequality gives P(theta0 outside C)<=.05. Conditional normalized likelihoods
apply this argument conditional on the held-out retained counts. This is the
established [split likelihood method](https://arxiv.org/abs/1912.11436).

The generating parameter belongs to all **54/54 matched sets**. With three
repeats per case this is not an empirical validation of 95% coverage. We have
not computed certified intervention-image widths of these sets, so point
accuracy must not be advertised as a useful prediction guarantee.

Misspecified-model membership is diagnostic only. With wrong perturbation
strength, base parameters belong to 6/6 normalized sets, 1/6 retention sets and
0/6 scheduled sets. Membership in a weak set is not proof of robustness; these
counts do not estimate coverage of a correctly specified parameter.

## Decision and next research question

Stop promoting the normalized design as a better estimator. Keep the
identifiability counterexamples and strong controls. The strongest remaining
question is: **what minimum calibration makes an intervention prediction
identifiable and precise when observation loss and intervention strength are
uncertain?**

A paper-oriented next phase must establish a substantial result about the
intervention target, rather than recover every nuisance parameter. First derive
identified sets when perturbation strength lies in a declared interval and
loading is uncertain; separate exact nonidentification from numerical looseness.
Then construct certified predictive bounds and compare standard information-based
design and retained-absence controls under an explicit cost model. Freeze a
fresh evaluation only after that mathematical claim survives prior-art review.
A biological claim additionally needs a measured held-out perturbation with a
validated observable. Neither that validation nor equation-level priority is
established here. This round supplies no new learning architecture.

A [post-screen exact-data deduction](UNKNOWN_PERTURBATION_THEORY.md) makes this
question more specific: unequal state susceptibility can identify unknown
perturbation strength from two absolute snapshots. Equal susceptibility instead
admits an exact 3.975-point target ambiguity in the constructed case. Both claims
are tested; useful finite-data calibration and originality remain unestablished.

## Reproduction

    python study_finite_loss.py
    python verify_finite_loss.py
    python -m unittest discover -s tests -p test_finite_loss.py -v

The verifier does not refit. Small tests compare numerical backends at domain
extremes, check state-dependent loss and exactly enumerate a binomial example of
the split-likelihood guarantee. All failed research gates remain in the record.
