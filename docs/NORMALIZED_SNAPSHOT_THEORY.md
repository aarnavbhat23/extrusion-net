# Candidate paper result: calibration of normalized kinetic snapshots

Status: algebraic development note with high-precision illustrative checks.
Priority over the identifiability literature is unestablished. No new biological
data were fitted. The result changes the observation contract and does not
contradict the earlier known-step-rate inverse.

## The observation contract is essential

Each independent episode starts from a known loading state. Observe at an
independent Exp(gamma) time with known gamma. Retain the spatial state only when
the motor is present, and record a fixed number of such retained states. The
data are a normalized, time-pooled snapshot histogram: total initiated episodes,
absences and individual observation timestamps are unavailable. Arm stepping
rates are unknown. Actual time-stamped trajectories contain more information
and are not covered by the impossibility statement.

## General normalized-resolvent ambiguity

Let T be a finite transient Markov subgenerator with nonnegative off-diagonal
rates and nonpositive row sums, and alpha the initial row distribution. The
unnormalized state probabilities at the independent observation clock are

    p_gamma = gamma * alpha * (gamma I - T)^(-1).

For c>=1 define

    T_c = c T - (c-1) gamma I.

This is another valid subgenerator, with transitions sped up by c and killing
rates changed from k_i to c*k_i+(c-1)*gamma. Since

    gamma I - T_c = c*(gamma I-T),

we have p_gamma,c=p_gamma/c. Thus normalization conditional on retention erases
the difference exactly. Any common fixed linear observation/emission map also
preserves this proportionality. Additional coarse-graining cannot recover the
lost information. This is elementary resolvent algebra and an application of
classical identifiability reasoning, not a new general Markov formalism.

## Six-parameter capture family

Extend the existing model to unknown positive speeds vL,vR, positive mu,h, and
capture probabilities pL,pR. Keep loading, barriers at three steps and geometry
known. The ambiguity remains inside this family:

    vL' = c*vL,  vR' = c*vR,
    mu' = c*mu+(c-1)*gamma,
    h'  = [c*mu*h+(c-1)*gamma]/mu',
    pL' = pL, pR' = pR.

The parameters are held fixed when applying a specified intervention that
multiplies both free and captured removal rates by f. Put R_f=vL+vR+f*mu+gamma.
The ratio of intervention snapshot probabilities P[1,0]/P[0,0] is vL/R_f,
regardless of normalization. Under the transformed parameter it becomes

    vL / [R_f + (f-1)*(c-1)*gamma/c].

For c>1 and f!=1 these differ. Consequently identical calibration histograms
can imply different intervention distributions even with unlimited calibration
samples. This is structural nonidentifiability under the stated data contract,
unlike the earlier finite-data confidence-set ambiguity.

## Statistical consequence

Choose any two such mechanisms whose scalar target values differ by Delta.
Their retained calibration samples have the same law at every sample size.
The triangle inequality implies that every estimator has maximum expected
absolute error at least Delta/2 across the pair. Any interval with coverage at
least 1-alpha at each mechanism must contain both target values with probability
at least 1-2*alpha under their common calibration law, hence have width at least
Delta on that event. These are standard indistinguishability/two-point arguments.
They cease to apply if timestamps, retention counts or calibrated speeds are
included in the data.

## Two clocks restore identification for this family

Let q_j[a,b] be the normalized present-state probabilities at distinct known
gamma_1,gamma_2. Define

    a_j = q_j[1,0]/q_j[0,0],  b_j = q_j[0,1]/q_j[0,0].

For the specified capture family,

    vL = (gamma_2-gamma_1)/(1/a_2-1/a_1),
    vR = vL*b_1/a_1,
    mu = vL/a_1-vL-vR-gamma_1,
    pL = 1-q_1[4,0]/(q_1[2,0]*a_1^2),
    pR = 1-q_1[0,4]/(q_1[0,2]*b_1^2).

For pL>0, let D=q_1[3,0]/q_1[2,0]-(1-pL)*a_1. Then

    h = (vL*pL/D-vR-gamma_1)/mu.

Use the transposed expression if only pR>0. If both captures are zero, h is
irrelevant to all removal interventions. These formulas prove global
identification on this family. They do not prove identification of an arbitrary
hidden generator from two emission distributions.

Alternatively, one calibrated speed or the absolute retention probability
breaks this particular ambiguity. Individual timestamps are another important
information source to test; a second acquisition must not be sold as necessary
when inexpensive existing metadata suffice.

## Exact recovery from one absolute-retention snapshot

If initiated episodes and absences are observed, P[0,0] is an absolute
probability. Set R=gamma/P[0,0], a=P[1,0]/P[0,0], b=P[0,1]/P[0,0]. Then
vL=R*a, vR=R*b and mu=R-vL-vR-gamma. The capture formulas above then recover
the remaining relevant parameters using one clock. A second clock is not
necessary when the denominator is observed accurately.

## Finite-sample speed diagnostic

For each clock, condition on the total number of observations in states (0,0)
and (1,0). The former count is binomial with success probability
w=q00/(q00+q10). Construct a Clopper-Pearson interval at error alpha/2 for each
clock and transform it monotonically to z=w/(1-w)=q00/q10. By a union bound,
both z values lie in their intervals with probability at least 1-alpha,
including after averaging over the conditioned totals. Zero totals give [0,infinity].

Subtract intervals to obtain [dL,dU] for z2-z1. Under the family the difference
is (gamma2-gamma1)/vL>0. If dU<=0, the confidence set is empty. Otherwise invert
its positive part: [(gamma2-gamma1)/dU, (gamma2-gamma1)/dL] if dL>0, or an
unbounded upper endpoint otherwise. This is a speed confidence set, not an
intervention prediction interval. Exact binomial coverage refers to the
mathematical construction; software evaluates quantiles numerically.

The guarantee conditions on retained observations. In the finite-data study,
initiated counts are available only to the resource ledger for the two-clock
strategies; their likelihood discards the count information. If that count is
actually known to the analyst, the full-retention comparator uses information
the other strategies deliberately discard. The impossibility theorem must not
be applied to that richer observation contract.

## Evidence and priority

Boundary tests independently assemble state generators and check the inverse,
the normalized ambiguity, intervention separation and one-clock absolute recovery.
The [frozen finite-data study](CLOCK_CALIBRATION_RESULTS.md) tests practical
conditioning and fails its continuation gate. Its formal speed intervals are
unbounded in all 48 two-clock fits. Structural identification is not a precision
or efficiency result.

Identifiability and experimental calibration are established research areas;
see [the kinetic identifiability framework](https://academic.oup.com/bioinformatics/article/35/24/5216/5514038),
and [Schlicke's snapshot time-scale analysis](https://arxiv.org/abs/2610.07341).
The latter uses stationary reaction-diffusion spectra, unlike the killed-chain
random-time snapshots here. That difference does not establish priority for
our elementary resolvent identity. A complete equation-level priority review
and external mathematical review remain required.
