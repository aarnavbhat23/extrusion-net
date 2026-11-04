# What a pulse-order signature cannot attribute

## 1. Two mechanisms with identical observations

Use row-state convention and initial state p0=(1,0). Let

    Q_A = [[-a,a],[0,0]], Q_B = [[0,0],[b,-b]], P_u=exp(tau Q_u).

In model M, state 1 means biological binding. In model O, state 1 means an
optical dark state; biological binding never occurs. Give both the same primary
emission distribution, e.g. Bernoulli detection probabilities (1,.4).

For any pulse word u1,...,um and binary observations y1,...,ym,

    L(y | u) = p0 P_u1 E_y1 ... P_um E_ym 1,

where E_y is the diagonal emission matrix. All factors are identical in M and
O, so the complete observed laws agree, not just their means. If a randomized
policy chooses each next pulse from the preceding observations, induction gives
the same transcript law in both worlds. Adaptive pulse selection, finer sampling
and more molecules cannot distinguish these model interpretations. Equal-prior
binary classification error is at least 1/2. A Bayesian preference can choose a
world but does not add information to the observation law.

This is an intentionally adversarial model class: the optical state is permitted
to respond to the same pulses. If an independently validated experimental
restriction rules that out, the counterexample may no longer apply. That
restriction cannot be inferred from an order effect alone.

## 2. One additional measurement, with its exact assumption

For unit pulses, x_AB=(1-exp(-a))*exp(-b), x_BA=1-exp(-a), and

    Delta=x_BA-x_AB=(1-exp(-a))*(1-exp(-b)).

Suppose an independently calibrated reporter has response f+kappa*x to
biological occupancy and no optical-state sensitivity. Its population contrast
d then satisfies d=kappa*Delta+u, with |u|<=2 epsilon if each arm's probability
offset is at most epsilon. A shared f cancels. If kappa is known and epsilon=0,
Delta=d/kappa; if kappa lies in [kL,kU], kL>0, invert over that interval and u.

With d known only to [dL,dU], the resulting contrast-level set is

    { (d-u)/k : d in [dL,dU], |u|<=2 epsilon, k in [kL,kU] } intersect [-1,1].

Its endpoints are the minimum/maximum of the four endpoint ratios. The bounds
are sharp for this reduced contrast model. They need not be sharp after adding
separate arm-level probabilities, calibration observations or kinetic constraints.
An empty interval signals incompatibility, not identification. If kappa=0 is
admissible, one cannot use this division formula to claim identification.

For n independent Bernoulli endpoints per arm, Hoeffding gives simultaneous
contrast coverage at least .95 with radius sqrt(log(40)/n). Inverting the
whole interval preserves coverage under the stated fixed calibration bounds.
Those bounds are assumed, not estimated by this screen; uncertain calibration
coverage would need its own error allocation. This is fixed-time inference.

## 3. Defeating the extra reporter

Give the optical state the same second-reporter emissions f+kappa*x as the
biological state, with the same joint emission law across both channels. The
full-law equivalence now extends to the added reporter. Even two channels do
not establish binding without justified selectivity/cross-talk constraints.

Likewise, an A-versus-baseline direct state reporter already separates M and O
under the ideal selectivity assumption. It measures existence of binding rather
than the order-effect Delta. A successful order test therefore cannot establish
that pulse ordering is necessary or efficient for this binary attribution task.
The earlier claim that order alone exposes biological memory must be narrowed.

## 4. Originality boundary

This elementary controlled-HMM construction is not a new equivalence theorem.
See the explicit [equation/protocol comparison](OPTICAL_MEMORY_PRIOR_ART.md).
The useful output is an auditable failure boundary and a precise missing assay
requirement. A substantive new paper would need an achievable reference
measurement or independently justified intervention restriction, and a useful
target-prediction guarantee beyond conventional identifiability analysis.
