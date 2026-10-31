# Nature-inspired direction: distinguish regulated turnover from observation loss

Date: 2026-10-29. The [structural protocol](LOSS_PERTURBATION_PROTOCOL.md) was
committed as `062bef1` before evaluation. This is an exact-probability feasibility
screen, not finite-sample success, a learned architecture, or experimental validation.

The biological inspiration is competing processes with different responses to
interventions. A motor can disappear because it is removed, while its observed
signal can also disappear for another reason. Treating both as one kinetic rate
can fit a baseline while predicting the wrong intervention. Conversely, adding
an unconstrained hidden process can make even ideal measurements insufficient.

## Model and observed results

Use parameters (vL,vR,mu,pL,pR,h,delta). With known relative speed s, known
turnover multiplier f and observation rate gamma, free loss is f*mu+delta,
captured loss is f*mu*h+delta, and arm speeds are s*vL,s*vR. Delta is a shared,
unaffected loss. Loading and capture structure remain known and invariant.
This model does not assert that any real cofactor or reagent changes only f.

| Design | Parameters varied in diagnostic | Local rank in all 3 mechanisms | Interpretation |
|---|---:|---:|---|
| Three relative speeds, one clock | 7 | 5 | Two ambiguities remain |
| Three turnover levels, delta known zero | 6 | 6 | Calibration can work under this stronger assumption |
| Three turnover levels, delta unknown | 7 | 6 | A scale/loss ambiguity remains |
| Two clocks, one turnover level | 7 | 6 | Regulated and unaffected losses remain confounded |
| Two clocks at baseline + one turnover perturbation | 7 | 7 | Explicit inverse recovers relevant parameters |
| Absolute retention at two turnover levels | 7 | 7 | Explicit inverse also recovers relevant parameters |

Ranks agree at both relative finite-difference steps (1e-4 and 5e-5).
Numerical rank is only a local diagnostic. The actual evidence for ambiguity is
an explicit observationally equivalent pair; evidence for global recovery in
the sufficient designs is the algebra below. Maximum reconstruction error over
six sufficient-design checks is **4.67e-15**. One-sided and no-capture boundary
tests are separate. When both capture probabilities vanish, h is irrelevant.

For each mechanism, a pair with identical baseline observations at both clocks
still differs under doubled regulated removal. Normalized target total variation
is **.01939, .01017 and .00818**. This ambiguity remains even if absolute retention
is measured at baseline: separating clocks does not by itself separate the two
biological/measurement causes of disappearance.

## Explicit indistinguishable pairs

At fixed gamma, scaling speeds and regulated removal by c, leaving capture
probabilities and h unchanged, and setting

    delta' = c*delta+(c-1)*gamma

gives T_u'=c*T_u-(c-1)*gamma*I for every relative speed and turnover input u.
The retained probabilities scale by 1/c. All their normalized histograms agree.
Thus arbitrarily many such steady-condition perturbations do not remove this
particular ambiguity when unaffected loss is unknown. Transition-rate ratios,
and consequently any constraints expressed solely through those ratios, survive
the same rescaling; thermodynamic consistency alone need not fix the scale.

Separately, at f=1, for feasible epsilon define

    mu' = mu-epsilon, delta' = delta+epsilon,
    h' = (mu*h-epsilon)/(mu-epsilon).

Both effective losses remain unchanged. Every baseline clock distribution,
including absolute retention, is identical, but responses at f!=1 differ.
The screen uses epsilon=.05 and c=1.6; all parameters remain valid. These are
restricted-family counterexamples, not evidence that arbitrary experiments fail.

## Sufficient three-condition inverse

Use normalized observations at (f,gamma)=(1,.1),(1,.3),(2,.3), with s=1.
The two baseline clocks identify vL and vR using the previously derived
[early-state ratio inverse](NORMALIZED_SNAPSHOT_THEORY.md). At any condition,
with a=q10/q00 and b=q01/q00,

    k0(f) = vL/a-vL-vR-gamma,
    pL = 1-q40/(q20*a^2), pR = 1-q04/(q02*b^2).

If pL>0, put D=q30/q20-(1-pL)*a; then

    kc(f) = vL*pL/D-vR-gamma.

Use the right-arm expression if only pR>0. Now

    mu = k0(2)-k0(1), delta = 2*k0(1)-k0(2),
    h = (kc(1)-delta)/mu.

These recover the relevant parameters globally on the stated positive-rate
family. For absolute probabilities at two turnover levels, obtain the speeds
from gamma/P00 and the early-state ratios at the first condition; the same
remaining equations apply. This retains the stronger absolute-retention control.

## General design insight, derived after the screen

Suppose the generator has the form

    T_u = Q_u(theta) - diag(sum_r phi_r(u)*kappa_ir),

where Q_u is conservative and scaling theta by c scales Q_u by c. A uniform-loss
rescaling can be constructed whenever a shared beta satisfies

    sum_r phi_r(u)*beta_r = gamma_u for every observed condition u.

Set kappa'_ir=c*kappa_ir+(c-1)*beta_r for every state i. Subject to parameter
feasibility, this gives the same normalized-resolvent equivalence. Conversely,
this shared coefficient-offset construction requires that linear system to be
solvable. Thus gamma lying in the column space of the loss-design matrix is an
audit for this particular ambiguity. Positive coefficients admit sufficiently
small feasible perturbations; boundary constraints require separate checking.

An unknown background-loss column of ones accommodates any constant gamma,
even with many biological interventions. Repeating the same biological condition
at different known gamma values removes this offset construction. It does NOT
prove full identification: the loss-decomposition pair above is a counterexample.
The three-condition inverse supplies the missing argument in our small family.
This design observation is not claimed as a new general identifiability theory.

## Direction and novelty boundary

The selected question is: **What is the smallest calibrated experiment set that
separates regulated kinetics from observation loss well enough to predict an
unseen perturbation?** The next phase should compare finite-sample prediction
certificates, timestamp likelihoods, and absolute-retention controls, with
misspecified turnover strength and dropout. Use fresh mechanisms and a locked
error/coverage criterion; do not treat this exact-data success as that test.

The nature-inspired model must earn value through transfer across environments.
It is not yet a new way to learn representations, and a seven-parameter inverse
alone is unlikely to justify an 8–9-level paper. A broader useful result or an
independently validated prediction is required before upgrading novelty.

Relevant primary precedents:

- [Shah and Tortora et al.](https://doi.org/10.1016/j.molcel.2026.03.037) already
  study cofactor dosage and extrusion rates; their
  [simulation code](https://github.com/Fudenberg-Research-Group/2026_shah-tortora_extrusion-rate)
  is a potential independent comparator. It does not establish our selective
  f multiplier or observation contract.
- [Finding and Breaking Lie Symmetries](https://www.mdpi.com/2073-8994/12/3/469)
  is an established biological-identifiability approach. Breaking a displayed
  symmetry alone does not prove identification; our inverse is essential.
- [Force-dependent extrusion and mechanical response](https://doi.org/10.1103/qfdv-vhv1)
  already includes hysteresis. Adding force feedback is not sufficient novelty.

Reproduce with `python study_loss_perturbations.py --verify` and the full test
suite. [All outcomes](../results/loss_perturbations.json) retain all 18 design/case
combinations. Probabilities are checked against separately assembled matrix and
flow implementations; saved-result replay is deterministic, not external review.
