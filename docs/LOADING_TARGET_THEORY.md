# Clock-specific prediction despite unknown loading and loss

2026-11-01. Restricted identification and interval-image results. Linear programs,
exact binomial intervals, resolvent algebra and partial identification are
established tools. The application question is which calibration is needed for
one intervention target, even when kinetic parameters remain unidentified.

## A loading gauge for complete snapshots

At a fixed exponential clock gamma, let loaded-state subgenerator be
T_f=Q-f*K-delta*I, where Q has nonnegative transition rates, K is nonnegative
diagonal regulated loss and delta is unaffected observation loss. Each initiated
episode loads with common probability ell and initial state distribution alpha.
The absolute observed present-state law is

    p_f = ell*gamma*alpha*(gamma*I-T_f)^(-1).

For s>=1 with s*ell<=1, choose Q'=s*Q, K'=s*K,
delta'=s*delta+(s-1)*gamma and ell'=s*ell. Then

    gamma*I-T'_f = s*(gamma*I-T_f), so p'_f=p_f.

This holds for every f at the SAME clock, including all common linear emissions
and the complementary absence probability. Thus absolute retention alone does
not identify loading/rate scale in this enlarged model. Yet every intervention
snapshot at this same clock is invariant along this particular gauge. Parameter
nonidentification does not imply nonidentification of every intervention target.
At a different observation clock the gauge generally fails; changing clocks can
therefore expose ambiguity that no number of same-clock f arms removes.

This extends the previous normalized-snapshot example by allowing failed loading
and separate optical loss. It does not prove that all observational equivalences
preserve interventions, or that a general unknown generator is identified.
Independent full-generator tests verify complete equality at three multipliers
and separation at a different clock for a concrete example.

## Scalar intervention target when susceptibility is equal

Assume h=1: all loaded states have total loss a=mu+delta at baseline and
a+(f-1)*mu under calibration intervention f>1. Let the target g=3 use baseline
loading and the same gamma=.3. Rates mu,delta are nonnegative; ell lies in
[ell_min,ell_max] with ell_max<=1. Let R0,R1 be absolute presence probabilities.
The new parameter class has nonnegative hazards without the earlier numerical
optimizer's finite rate box. Sharp endpoint mechanisms need not belong to that
older fitting box; imposing it would require extra constraints and a new image.
Under common loading,

    U=1/R0=(gamma+a)/(ell*gamma),
    V=1/R1=(gamma+a+(f-1)*mu)/(ell*gamma),
    1/Rg=U+(g-1)/(f-1)*(V-U).

The unknown loading factor cancels. At known f this identifies the target
presence from exact probabilities despite nonidentified mu,delta,ell. For
f in [fL,fU], put tL=(g-1)/(fU-1), tU=(g-1)/(fL-1), introduce W=(g-1)*mu/(ell*gamma),
and impose

    V>=U, tL*(V-U)<=W<=tU*(V-U),
    U>=1/ell_max, W<=(g-1)*(U-1/ell_max).

Probability intervals constrain U and V to reciprocal intervals. Minimize and
maximize Z=U+W over this polytope. Since target absence is 1-1/Z, the extrema
give its exact identified interval for this probability-rectangle model.

Necessity follows by substitution and delta>=0. For sufficiency, take any
feasible point and choose ell=ell_max,

    mu=ell_max*gamma*W/(g-1),
    delta=gamma*(ell_max*U-1)-mu,
    f=1+(g-1)*(V-U)/W when W>0.

These are physical parameters, reproduce the two observed probabilities and
attain the stated target. W=0 implies V=U and mu=0, with f irrelevant. The closed
nonnegative-rate class includes that boundary; with strict positivity it may
instead describe an attained limit. The lower loading bound does not tighten the
image because the construction uses the permitted maximum loading.

## Bounded loading drift: a sharp four-variable construction

Now calibration-arm loading is kappa*ell while target loading remains ell.
Assume kappa in [kL,kU] subset (0,1]. Set J=kappa*V. The exact constraints become

    kL*V<=J<=kU*V, J>=U,
    tL*(J-U)<=W<=tU*(J-U),
    U>=1/ell_max, W<=(g-1)*(U-1/ell_max).

Again optimize U+W. Every physical mechanism maps into the polytope. Conversely,
take ell=ell_max, recover mu/delta as above, kappa=J/V and
f=1+(g-1)*(J-U)/W. Since kappa<=1, treated loading is also <=1. This proves
sharpness and yields explicit endpoint mechanisms. The construction needs extra
feasibility treatment if loading can increase above baseline (kappa>1). It does
not model arbitrary state-specific loading or observation errors.

## Finite-sample confidence and certified computation

Two independent binomial retained/initiated counts give standard exact probability
intervals, each at error .025. A union bound covers both true probabilities with
at least .95 probability under the model. Mapping their rectangle through the
entire feasible polytope therefore gives target coverage at least .95. Sharpness
of that image does not mean the confidence rectangle is the shortest possible
confidence procedure. Quantiles and LPs are evaluated numerically.

Every saved optimum has a feasible primal point and a feasible dual point with
matching objectives, checked separately without running an optimizer. This
certifies global extrema of the declared finite LP within numerical tolerance.
Zero lower probability bounds use conservative [0,1], not reciprocal truncation.
Solver-incompatible confidence sets also return [0,1] with an explicit status.
Their width is not an exact-data identified width; incompatible exact-data sets
are marked empty rather than presented as a one-unit identified interval.

Four scheduled-time binomial intervals provide a conventional alternative.
With L=log(ell), A=a/gamma, D=(f-1)*mu/gamma, W=(g-1)*mu/gamma, log probabilities
at times 1/gamma and 2/gamma give linear constraints. Target presence is
exp(L)/(1+A+W). Separate extrema give an outer interval; 16 fixed L slices
reduce its looseness with a recorded one-sided slack bound. The timestamp
interval remains an outer approximation and uses different observation durations.
This comparison cannot prove superiority over optimal timestamp inference.

## Failure boundaries and priority

Condition-varying loading invalidates the common-loading interval, as demonstrated
by its stress failure. The ratio-robust interval requires a valid kappa contract;
the narrow interval in the study is externally supplied, not estimated from the
same survival counts. Unknown target loading would require another uncertainty
variable. State-dependent susceptibility h!=1 invalidates the scalar reciprocal
formula; a separate full-generator test demonstrates that failure. The general
loading gauge above and the restricted scalar target have different scopes.

Established partial-identification methods already optimize constraints under
measurement error; see [Finkelstein et al., 2021](https://proceedings.mlr.press/v161/finkelstein21b.html)
and [automated sharp bounding](https://pmc.ncbi.nlm.nih.gov/articles/PMC11566246/).
They prevent presenting LP sharpness or uncertainty projection as a new general
framework. The possible narrow contribution is this clock-specific kinetic
target, its loading-drift reduction and the resulting calibration requirements.
Priority for that application result remains unestablished and requires a
complete comparison and independent mathematical feedback before submission.
