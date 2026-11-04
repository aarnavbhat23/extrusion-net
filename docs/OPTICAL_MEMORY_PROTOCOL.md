# Biological versus optical memory: frozen comparison

Date: 2026-10-30 (local). This follows, rather than replaces, the failed
2026-10-21 pulse-order screen. The new target is attribution of an order effect,
not discovery of noncommutativity or another schedule family. No measured
biological data are used to select or evaluate this screen.

## Smallest counterexample and proposed measurement

Two states, initially state 0. A makes 0 -> 1 at rate a; B makes 1 -> 0 at
rate b. Each pulse lasts one time unit. In a biological model state 1 is bound;
in an optical model state 1 is dark and the molecule is never bound. Give both
models the same primary emission probabilities (1, .4), transition matrices
and initial state. Their complete primary observation laws agree under every
input sequence, including adaptively selected sequences. Biological occupancy
differs. This is a mechanism-label ambiguity, not a failure of numerical fitting.

An additional calibrated reporter has detection probability f+kappa*x for
biological bound occupancy x, with f=.1 and kappa=.7. In the optical-only model
it has probability f. The target is Delta=x_BA-x_AB, not a gene-expression or
general-learning claim. Permitted arm-specific probability offsets have absolute
value <=epsilon=.01. Calibration supplies kappa in [.6,.8] and excludes optical
cross-talk. The effect of relaxing that last assumption is explicitly tested.

## Larger fixed screen, not yet evaluated when frozen

Six mechanisms: a in [.2,.6,1.5], b in [.3,1.2], pulse duration 1. For each,
sample n=500 and n=2000 independent endpoints per arm, 2,000 independent
experiment replicates per setting, seed 2026103001. Equal numbers of endpoints
do not imply equal instrument costs across reporters.

Record primary-only order tests for both mechanisms; calibrated second-reporter
order tests; and a simpler A-versus-baseline reporter control with the same
2n endpoints. The simpler control tests existence of binding, not the same
order-effect estimand; it prevents treating pulse order as necessary to
distinguish these particular two models.

Use a fixed-time Hoeffding difference radius sqrt(log(2/.05)/n). Expand by
2*epsilon and invert over the entire calibration interval to bound Delta.
Under biological alternatives use offsets that attenuate the contrast; under
optical nulls use offsets that maximize a positive contrast. Also simulate an
out-of-contract optical reporter whose state sensitivity equals kappa: it can
reproduce the biological reporter law and must not be called a validated null.
Report Wilson intervals for rejection/coverage/power estimates.

Technical gate: all schedule-law counterexamples agree to 1e-12; minimum
in-contract coverage >=.94 in every setting; optical-null rejection Wilson
upper bound <=.075; at n=2000, >=.80 order-test power in at least four of six
mechanisms. Failures remain failures. Even a pass is only conditional method
feasibility, not novelty, biological validation or an advantage over the simpler
control. Never retune after seeing these results.

## Identification and failure boundaries to verify

Delta=(1-exp(-a))*(1-exp(-b)). The reporter contrast is
d=kappa*Delta+u, |u|<=2*epsilon. Nonzero kappa and a bound on u are essential.
If kappa can be zero or the optical state can produce the same emissions in
both channels, increasing n or changing the pulse sequence cannot identify
biological attribution. Compare analytic occupancy with matrix exponentials;
enumerate short observation strings to verify full-law equivalence.

Originality audit must compare HMM observational equivalence, fluorescence
blinking, observability/experimental design and the repository's prior pulse
work. No claim that the algebra, confidence bound or second-reporter principle
is new. A substantive paper requires an assay-specific justified restriction
on optical cross-talk plus independent state measurements.
