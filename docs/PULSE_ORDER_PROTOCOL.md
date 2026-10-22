# Endpoint-only pulse-order screen — frozen 2026-10-21

Question: can equal-exposure schedules detect a hidden protective state while
remaining insensitive to arbitrary static molecular frailty and low-order
force-by-age hazard? This is a proposed experiment/model-falsification tool,
not a new learned architecture or measured CTCF result. The evaluated CTCF
workbook will not be reused to select this candidate.

## Null and design

For a molecule z, let its one-event hazard under condition u at age t be
lambda_z(u,t)=sum_(j=0)^d c_(z,u,j)*t^j, nonnegative over the experiment. Draw z
from the same population independently in both randomized arms. Integrating
the hazard gives endpoint survival. Arbitrary dependence/heterogeneity among
the molecule-specific coefficients is allowed. The schedule cannot change z.

Use a binary Thue–Morse word of length 2^(d+1) and its complement, with equal
slot durations. Prouhet moment cancellation gives equal integrated hazards
for every z up to degree d, hence equal mixture survival. AB/BA cancels constant
hazards, ABBA/BAAB also cancels linear ageing, and the eight-slot pair also
cancels quadratic ageing. This is classical moment balancing, not a newly
invented sequence. Equal total duration T=4 and equal A/B exposure across all
schedule lengths; report switch counts rather than claiming equal hardware cost.

## Fixed simulator screen

Use row-vector two-state transient generators, initial distribution (1,0).
State 0 is unprotected, state 1 protected. Under A: capture k, release .05,
event hazards .2 and .2*h. Under B: capture .05, release k, event hazards 1.2
and 1.2*h. Cases k=.15/.5/1.5 by h=.1/.4 (six total), not fitted to biological
data. Propagate exactly by matrix exponentials; verify against a separate
three-state event simulation that includes the absorbing outcome.

Compare the two-, four- and eight-slot complementary schedules at fixed T=4.
Draw 1,000 independent binary endpoints per arm, 2,000 experiment replicates,
seed 2101. Primary four-slot family. Known systematic endpoint corruption budget
r=5 per arm: a Hoeffding difference test rejects if the observed absolute
frequency difference exceeds
sqrt(.5*(1/nA+1/nB)*log(2/.05)) + r/nA + r/nB.
Simulate worst-direction corruption: away from zero under null, toward zero
under alternatives. This is a supplied corruption bound, not an estimated
assay error rate. Also report the ideal two-sided Fisher test with no corruption.
Each schedule family is a separate diagnostic; no family-selected combined
claim at level .05. A future multi-family test requires multiplicity control.

Null screen: 100 fixed-seed random finite frailty mixtures with 32 components,
coefficient arrays shape (32,2,2), c0~Uniform(.01,1), c1~Uniform(0,.2),
Dirichlet(1) weights. Three sequence families, same endpoint sampling and
replicates. Four-slot success: all 100 exact survival contrasts below 1e-12,
maximum simulated robust null rejection <=.065, and >=.80 robust power in at
least four of six hidden-state cases. Failures stop an experimental-efficiency
claim. Report a force-by-quadratic-age counterexample beyond the four-slot null,
and an unobservable hidden-state example with equal event hazards across states.
Nonzero contrast rejects the specified scalar-hazard null; it does not uniquely
identify CTCF capture or establish a particular number of hidden states.

## Novelty boundary

Force-dependent CTCF barriers: Davidson et al., https://doi.org/10.1038/s41586-023-05961-5 .
Static/dynamic disorder in force spectroscopy already has theory and data:
https://pmc.ncbi.nlm.nih.gov/articles/PMC4941464/ . Drift-eliminating ABBA designs
are standard: https://itl.nist.gov/div898/handbook/mpc/section3/mpc3312.htm and
Swanson/Schlamminger, https://arxiv.org/abs/1009.1894 . Thue–Morse driving and
higher-order cancellation have extensive prior art, e.g.
https://arxiv.org/abs/2101.07065 . The candidate is the restricted endpoint-survival
assay with an explicit frailty/age null, bounded observation corruption and an
observable-signal/power tradeoff. Search results do not establish priority.
