# Loading-drift follow-up — frozen after the common-loading failure

2026-11-01. The first common-loading screen passes its matched gate, but treated
loading .8*baseline produces 0% target coverage in two of three stress mechanisms.
Retain that failure. This follow-up changes the uncertainty contract explicitly;
it is not a reinterpretation of the first result or independent biological data.

Let baseline loading ell in [.6,1], treated loading kappa*ell, with kappa in [.8,1].
The target uses baseline loading ell. Keep h=1, gamma=.3, target multiplier g=3,
unknown calibration multiplier f in [1.8,2.2]. With U=1/R0,V=1/R1, introduce
J=kappa*V and W=(g-1)*mu/(ell*gamma). The probability-interval image is a polytope:

    kappa_L*V <= J <= kappa_U*V, J>=U,
    tL*(J-U) <= W <= tU*(J-U),
    W <= (g-1)*(U-1/ell_max).

Min/max U+W gives sharp target-absence endpoints. Construct every endpoint's
physical mechanism at ell=ell_max; kappa=J/V<=1 keeps treated loading <=1.
This proof needs a different feasibility treatment if kappa can exceed one.

Fresh mechanisms (mu,delta,ell)=(.2,.2,.75),(.4,.1,.85),(.65,.12,.95), true f=2,
true kappa=[1,.9,.8], N=24000 split across baseline/treated; 100 repeats per
mechanism/kappa; seed2026110103. Use exact binomial probability rectangles with
union error .05. Every method receives the same saved counts.

Compare (a) the previous common-loading assumption, (b) broad kappa [.8,1],
(c) an externally calibrated narrow kappa interval intersecting [kappa-.02,
kappa+.02] with [.8,1], and (d) exact known kappa. Narrow/exact controls have
extra information; their reference acquisition cost is unmeasured. No claim
that the narrow interval was estimated from these same survival counts.

Primary gates unchanged in form: per-method coverage Wilson lower>=.90 in every
mechanism/kappa cell, and median width<=.10 in >=2/3 mechanisms at EACH kappa.
Keep broad loading failure if it occurs. Report exact-data widths and endpoint
witnesses as structural uncertainty distinct from statistical noise. Mark empty
or zero-count cases and retain them. Independent primal/dual verification is
required; no local optimizer is a global certificate. The one-field spatial
reference experiment does not validate this kinetic loading contract.
