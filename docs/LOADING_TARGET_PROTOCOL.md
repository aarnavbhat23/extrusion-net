# Loading-uncertain intervention target — frozen 2026-11-01

Scope: the equal-susceptibility h=1 submodel of the existing capture process.
Independent failed loading occurs with common probability 1-ell in every arm;
loaded episodes start at the previously specified state. Baseline hazard a=mu+delta,
treated hazard a+(f-1)*mu, target hazard a+(g-1)*mu, gamma=.3, g=3. The same
ell lies in [.6,1] and unknown f lies in [1.8,2.2]. Speeds/capture do not affect
total survival when h=1. Do not generalize to h!=1 or state-dependent optical loss.

Theory to check: if R0,R1 are absolute presence probabilities at gamma, put
U=1/R0, V=1/R1, D=V-U. Target presence satisfies

    1/Rg = U+(g-1)/(f-1)*(V-U).

Common loading cancels from this prediction despite nonidentified biological
and optical hazards. Physical feasibility still depends on ell_max. Under
finite simultaneous probability intervals, introduce W=(g-1)*mu/(ell*gamma),
so tL*(V-U)<=W<=tU*(V-U), W<=(g-1)*(U-1/ell_max), V>=U.
Min/max U+W over this polytope gives the exact image of the interval-constrained
model (proof and endpoint witnesses required). Use the standard two-sided exact
binomial intervals with total union error .05. Zero lower probability bounds
give a conservative [0,1], never an arbitrary finite reciprocal cutoff.

Finite screen: three fresh (mu,delta,ell) mechanisms (.25,.1,.7),(.5,.15,.85),
(.7,.05,.95); f=2. Counts per design N=[6000,24000], split equally across baseline
and treated conditions, 100 independent replicates, seed2026110101. Absolute
snapshot observations retain failed loading and absence. Preserve all counts.
Compare uncertain-f intervals with the same estimator supplied accurate f=2
(an explicitly extra-calibration control), and a conventional timestamp design.
Both broad and accurately calibrated f controls use the same counts.

Timestamp control: separate initiated episodes, equal allocation across baseline/
treated and scheduled t=1/gamma,2/gamma (four groups total). Presence mean is
ell*exp(-r*t). Use four exact binomial intervals with union error .05. In variables
L=log(ell), A=a/gamma, D=(f-1)*mu/gamma, W=(g-1)*mu/gamma, the log-probability
constraints are linear, as are tL*D<=W<=tU*D and W<=(g-1)*A. Bound target presence
exp(L)/(1+A+W) conservatively by separate LP extrema of L and A+W. Report loss
from this relaxation; do not present it as a sharp timestamp optimum or prove
the proposed design superior to all timestamp inference. Account for durations.

Primary gates per design: minimum pointwise 95% coverage Wilson lower>=.90 in
all three N=24000 mechanisms; median target-absence interval width<=.10 in at
least two of three. A .10 width is a declared computational pilot tolerance,
not a biological effect threshold. Also report exact-data identified widths and
every failed gate. No repeated seed selection or changed tolerance is allowed.

Out-of-contract stress: N=24000, 100 replicates per mechanism, treated loading
.8*ell while fitting common loading. Keep matched-model coverage and stressed
diagnostic coverage separate. Compare to an intentionally fixed-ell=1 plug-in
point predictor on the SAME random-clock counts as a bias diagnostic, not as a
strong uncertainty method. Independent replay must check counts, interval images,
LP primal/dual certificates and source-generator survival probabilities.
