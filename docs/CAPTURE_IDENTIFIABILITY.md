# One-clock structural identifiability of the capture family

**Within the specified four-parameter family, one exact anchor snapshot
distribution identifies every parameter whenever capture is possible.** With no
capture, the unidentifiable parameter has no effect on any removal intervention.
The numerical ambiguity in finite samples therefore need not be a structural
alias between mechanisms. This is a small-family algebraic result, not evidence
of a novel biological mechanism or a general learning architecture.

## Assumptions and inverse

Keep the existing model: one motor loaded at (0,0), independent unit-rate arm
steps, known barriers three steps from loading, permanent capture on arrival
with probabilities pL,pR, removal rate mu before capture and mu*h afterward.
Observe at an independent Exp(gamma) time with known gamma>0. The 9x9 present
states and absent category are all distinguished. Assume mu>0, h>0, and exact
population probabilities P[a,b]. Loading, observation bias, rates and geometry
are known. These assumptions do not describe an ordinary uncalibrated contact map.

Write r=2+mu+gamma. Along the axis b=0 there is only the path in which the
right arm never steps. Competing-clock probabilities give

    P[0,0] = gamma/r
    P[2,0] = gamma/r^3
    P[4,0] = gamma*(1-pL)/r^5
    P[3,0] = gamma/r^3 * ((1-pL)/r + pL/(1+mu*h+gamma)).

The corresponding identities hold with axes exchanged. Hence

    r  = gamma/P[0,0]
    mu = r-2-gamma
    pL = 1-r^2*P[4,0]/P[2,0]
    pR = 1-r^2*P[0,4]/P[0,2].

If pL>0, put D=P[3,0]/P[2,0]-(1-pL)/r. Then

    h = (pL/D-1-gamma)/mu.

If only pR>0, use the transposed expression. The denominators are positive on
this model image. These formulas prove injectivity; no numerical optimizer,
Jacobian-rank assumption, or local-only argument is required. If pL=pR=0, h
never enters the process, so different h values give identical calibration and
intervention distributions. Every target at a specified removal multiplier is
therefore structurally identifiable even in that exceptional case.

## Why this does not settle finite-data prediction

The key beyond-barrier event P[4,0] decays as r^-5. Small empirical counts and
subtraction in the inverse amplify noise. Ratios of empirical frequencies can
even leave the parameter domain. The formulas are not a recommended estimator
or a finite-sample confidence interval. More observation clocks are not required
for structural identification; they may help precision and model checking.

The [joint likelihood audit](JOINT_LIKELIHOOD_RESULTS.md) retains pairs with
different intervention predictions after statistically valid joint testing.
Its split sample costs power, and its results are conditional on this family.
Neither that result nor this inverse establishes a sample-size lower bound for
every possible procedure. The next question is which experimentally accessible
measurement gives enough information to distinguish the relevant predictions.

The executable audit replays the identities with a separately implemented
70-digit forward flow for boundary cases, asymmetric capture, h below/above
one, and all six existing mechanisms. It reports expected and observed counts
for the informative states on the existing data. This is development analysis,
not a prospectively held-out biological validation or a priority claim over the
kinetic-identifiability literature. Modern measured cofactor models are a
required comparator: [Shah and Tortora et al.](https://doi.org/10.1016/j.molcel.2026.03.037).
