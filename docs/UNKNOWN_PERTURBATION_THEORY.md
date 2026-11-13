# State-dependent susceptibility as an internal calibration

2026-10-29. **Post-screen structural deduction, not finite-data confirmation.**
The wrong-dose stress in [the finite-loss study](FINITE_LOSS_RESULTS.md) motivates
estimating perturbation strength rather than fixing an inaccurate value. The
following elementary inverse works only under the specified shared kinetic
model. Priority is unestablished; it is not presented as a new general method.

Let the free and captured loss rates at baseline be a=mu+delta and b=mu*h+delta.
Under an unknown common multiplier f of regulated removal they become
c=f*mu+delta and d=f*mu*h+delta. Unaffected observation loss delta and state
susceptibility h are invariant across conditions. Assume known loading, absolute
snapshot probabilities at a known exponential clock, unchanged positive arm
speeds, and positive capture in at least one arm. The earlier local inverse
recovers speeds/capture and a,b,c,d from these exact snapshot distributions.

If mu>0, f!=1 and h!=1, then

    h = (d-b)/(c-a),
    mu = (a-b)/(1-h),
    delta = a-mu,
    f = 1+(c-a)/mu.

Substitution proves global recovery in this restricted domain. Thus unequal
state susceptibility can serve as an internal control for unknown perturbation
strength. Merely calling f uncertain does not always imply nonidentification.
The implementation is an exact-probability oracle, not a proposed noisy plug-in
estimator. Tests recover 24 combinations (including one-sided capture and both
increased/decreased turnover) from the independently assembled full generator.

The inverse explicitly divides by c-a and 1-h. Near no intervention or equal
susceptibility it becomes unstable; rare captured-state observations are another
finite-data bottleneck. Dropout that changes with state or intervention breaks
the shared-delta assumption. Unknown loading is not resolved by this result.

## Exact boundary ambiguity and a target interval

When h=1, a=b and c=d. Even infinitely many exact baseline and treated snapshots
at different observation clocks cannot split mu, delta and unknown f. Let
A=c-a>0, and suppose f lies in [fL,fU] with 1<fL<=fU. Feasible parameters obey

    mu=A/(f-1), delta=a-mu >= 0.

Intersect the multiplier interval with this physical constraint before using it.
For a new *known* multiplier g>1 the total hazard is

    r_g = a+(g-1)A/(f-1).

It is decreasing in f. Absence at clock gamma is r_g/(gamma+r_g), because h=1
makes removal state-independent. Therefore evaluating the two feasible interval
endpoints gives a **sharp exact-data identified interval** for that target.
Both endpoints are achievable, not just an outer relaxation. This interval
describes structural uncertainty, not a finite-sample confidence interval.

Example: a=.44, A=.256, f in [1.8,2.2], g=3, gamma=.3. The mechanisms
(mu=.32,delta=.12,f=1.8) and (mu=.213333...,delta=.226666...,f=2.2)
give identical complete baseline and treated snapshot distributions, but target
absence is .782609 versus .742857: a **3.975-percentage-point** unavoidable
exact-data span under this contract. The unit test checks equality at three
clock rates and the target separation with the full generator.

## Next decisive gate

Before another broad model search, ask whether this state contrast can support
useful predictive intervals with finite data near h=1. Compare a model that
estimates f with standard likelihood nuisance profiling and a genuinely measured
dose/retention control. Use a fresh frozen mechanism grid and a cost for dose
calibration. Require both coverage and width, and keep equal-susceptibility cases
as mandatory failures when the target interval exceeds tolerance. A certified
interval method must enclose the whole feasible parameter region; optimizer
endpoints alone do not establish this. This is a concrete next research question,
not a completed improvement in paper readiness or novelty.
