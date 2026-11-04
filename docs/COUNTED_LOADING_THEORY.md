# A counted reference is informative only with its detection contract

2026-11-01. Companion to LOADING_TARGET_THEORY.md and the frozen counted-reference
experiment. The probability intervals, ratio projection and union bound are
standard methods; the present application is a finite-information kinetic target.

Let loading be ell in the baseline arm and kappa*ell in the calibration arm.
A separate time-zero loading reference is positive with probabilities

    B0=q0*ell, B1=q1*kappa*ell,

provided its false-positive probability is zero and it measures loading before
loss. Independent reference episodes with known initiated denominators then
identify kappa=B1/B0 when q1=q0, even when their common sensitivity is unknown.
For q1/q0 in [qL,qU], only the interval

    kappa in [(B1/B0)/qU, (B1/B0)/qL] intersect (0,1]

is justified. A nonzero false-positive rate invalidates this multiplicative
relation unless it is explicitly accounted for.

Four exact binomial probability intervals, each error .05/4, jointly contain
the two retention and two reference probabilities with probability >=.95.
Their ratio image therefore contains the true kappa on that event. The existing
drift polytope contains the true kinetic point on the same event; its target
image covers with at least .95 probability. Independence is needed for the
declared binomial sampling within each arm, not for the union bound. Clustered
events from a specimen are not independent reference episodes by default.

The projection loses correlations from the joint four-probability model.
Although the retention target polytope is sharp for its supplied rectangle and
ratio contract, this composed confidence procedure is only a conservative outer
image. It is not a new exact joint likelihood method or optimal allocation rule.

## An exact counterexample to uncalibrated reporter sensitivity

Take gamma=3/10, f=2, target g=3 and ell=3/4 in both worlds. Consider:

| World | mu | delta | kappa | q0 | q1 |
|---|---:|---:|---:|---:|---:|
| A | 1/5 | 1/5 | 4/5 | 9/10 | 9/10 |
| B | 3/10 | 1/10 | 8/9 | 9/10 | 81/100 |

Both worlds give exactly the same four observed probabilities:

    R0=9/28, R1=1/5, B0=27/40, B1=27/50.

Their target absence probabilities differ: 35/44 versus 43/52, a gap 9/286
(about 3.15 percentage points). Thus arbitrarily many counts of all FOUR
observables do not identify the target if the sensitivity ratio is uncontrolled.
This is verified with exact rational arithmetic in the regression suite.
The observation contract, not the existence of a second reporter, makes the
reference informative. An independently calibrated sensitivity ratio or a
different reference operator is needed to break this equivalence.

## Priority boundary: the scalar identity is not a new kinetic law

Under common loading the exact scalar curve is

    R(f)=ell*gamma/(gamma+delta+f*mu).

Consequently

    R(1)/R(f)=1+(f-1)*mu/(gamma+delta+mu).

This is the same reciprocal-linear functional form as the classical
[Stern-Volmer relationship](https://old.goldbook.iupac.org/plain/S06004-plain.html)
(IUPAC S06004), with f-1 playing the role of the controlled quenching variable.
This is an algebraic equivalence of response curves, not a claim that the
physical experiments are identical. Two-point reciprocal interpolation and
loading cancellation cannot be presented as a new general mechanism.

[Finkelstein et al. (2021)](https://proceedings.mlr.press/v161/finkelstein21b.html)
already gives LP partial identification under discrete measurement error.
[Reisser et al. (2020)](https://doi.org/10.1038/s41598-020-58634-y) already separates
photobleaching and binding through illumination schedules. [Ciocanel et al.
(2024)](https://arxiv.org/abs/2307.15857) treats identifiable parameter combinations
for fluorescence models. The plausible contribution is a narrow, calibrated
target-prediction protocol and its sharp failure boundaries, subject to independent
priority review. The current experiments do not establish a new learning paradigm.
