# Aggregate survival identity: development refinement after interval failure

All six continuous-domain runs enclosed their generating mechanisms but missed
.02 precision. The largest target category width was the absent category in
every run. Original outcomes and the failed precision gate remain unchanged.
This amendment is a theory-driven reanalysis of those same data, not independent
confirmation and not a new action policy.

Derive survival to Exp(gamma) from the first successful capture time S. Before
S the removal rate is mu; after S it is mu*h. Write s=mu+gamma and
K=E[exp(-s*S); S finite]. Then

    P(present) = gamma/s * (1-K) + gamma/(mu*h+gamma) * K.

With two independent unit-rate arms and barriers three steps away,

    Ksingle = (1+s)^(-3)
    Kboth = 2/(s+2)^3 + 6/(s+2)^4 + 12/(s+2)^5
    K = [pL*(1-pR)+pR*(1-pL)]*Ksingle + pL*pR*Kboth.

This is a classical competing-risk/Laplace-transform calculation. Equivalently,
K=pL*H(3,0;s)+pR*H(0,3;s)-pL*pR*H(3,3;s), connecting the earlier reach kernel
to the randomized snapshot observable. No new biological mechanism is asserted.

By coupling, survival decreases with mu and h. Increasing a capture probability
increases survival when h<=1 and decreases it when h>=1. Therefore each rectangle's
exact absent-probability extrema occur at specified parameter corners, including
boxes spanning h=1. Evaluate those corners with outward arithmetic and intersect
the result with the existing interval-flow absent bounds. Keep every other
algorithm setting, input count and resource limit fixed.

Before applying the refinement, test the identity against the independent flow
solve and sampled-box inclusion across h below/above one. Then rerun enclosure
on the six saved calibration datasets, preserving original results separately.
Report containment, precision, compute and all box counts. Even if precision
improves, fresh mechanisms and a measured observation model remain required.
