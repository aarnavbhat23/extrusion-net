# What normalized fluorescence can identify

This is a restricted measurement-model deduction, not a new general
identifiability theorem or a validated description of this assay.

Let D(t) denote the difference between unbleached and bleached ROI signals and
C(t) the background-corrected reference. A reaction-limited decomposition is

    D(t)/D(0) = B(t) S(t),   C(t)/C(0) = B(t) V(t).

B is a shared optical attenuation factor, S the biological contrast of interest,
and V possible biological change in the reference. Source normalization observes
Y=S/V. Calling V=1 a reference assumption does not establish it experimentally.
Diffusion, spatial exchange and state mixtures may further change this operator;
the rate fitted to a curve is not automatically a molecular removal rate.

For any positive H with H(0)=1, the transformation

    (B,S,V) -> (B/H, S*H, V*H)

preserves BOTH raw signals, not just their ratio. Therefore weighting a raw
contrast objective cannot resolve this ambiguity on its own. The raw-contrast
loss is also exactly Z^2-weighted normalized-curve squared error, where Z=C/C(0).
This equivalence rules out describing that objective as a new learning mechanism.

## An explicit intervention ambiguity

In the exponential submodel B=exp(-b t), S=exp(-k t), V=exp(-v t), the raw data
identify lambda_D=b+k and lambda_C=b+v. With nonnegative rates, b can take any
value in [0,min(lambda_D,lambda_C)]. The shift

    (b,k,v) -> (b+z,k-z,v-z)

preserves complete noise-free D and C trajectories whenever rates stay feasible.
If a NEW known intervention multiplies only k by g>1, the normalized target is

    Y_g(t)=exp(-(g*lambda_D-lambda_C-(g-1)*b)*t).

It increases monotonically with b. The feasible endpoint values of b therefore
give a **sharp exact-data identified interval within this restricted submodel**.
Every endpoint corresponds to a feasible mechanism. This is structural
uncertainty, not a finite-sample confidence interval for the measured experiments.

Example, rates per minute: (b,k,v)=(.004,.030,.008) and (.009,.025,.003).
Both give lambda_D=.034, lambda_C=.012 and identical measured curves at every
time. Their predictions after g=2 differ at t=20 minutes. Allowing the complete
feasible b interval gives Y_g(20) in [exp(-1.12),exp(-.88)], a span of about
8.85 percentage points despite perfect baseline measurements.

An independent observation fixing b or v can break this particular ambiguity.
Merely collecting more baseline frames, fitting two exponentials, or reducing
optimizer error cannot. A biological claim requires a justified observation
operator and intervention selectivity in addition to good curve prediction.

This deduction narrows the next experiment: independently constrain reference
dynamics/attenuation, then test a perturbation whose biological target is measured.
It does not prove that this exact confounding occurs in the source experiments.
General FRAP reaction/diffusion identifiability has substantial prior work:
[Sprague et al.](https://doi.org/10.1529/biophysj.103.026765) and
[Ciocanel et al.](https://arxiv.org/abs/2307.15857). Priority of a useful extension
beyond those established problems remains unresolved.
