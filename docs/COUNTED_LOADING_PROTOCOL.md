# Counted reference calibration for the intervention target

2026-11-01. Freeze before outcomes. The earlier narrow loading-ratio interval
was supplied externally. Replace it with observed reference counts and account
for their finite-sample uncertainty. This is a simulated assay-design study,
not experimental validation of cohesin loading.

Keep the h=1 target and the three fresh mechanisms from LOADING_DRIFT_PROTOCOL:
(mu,delta,ell)=(.2,.2,.75),(.4,.1,.85),(.65,.12,.95); gamma=.3, calibration
factor f=2 in [1.8,2.2], target g=3. Treatment loading ratio kappa in {.8,.9,1}.
Each arm has 24000 initiated retention episodes. Independently count positive
loading-reference events among m initiated reference episodes in each arm,
m in {500,2000,8000,32000}. Reference detection sensitivity q0=q1=.9, with no
false positives; reference positive probabilities are q0*ell and q1*kappa*ell.
There are 100 repetitions/cell, seed 2026110105. Reuse retention counts across
reference budgets within each repetition; draw reference counts independently.

Construct exact binomial intervals on all FOUR probabilities, each error .0125,
so their joint rectangle covers with probability >=.95 by a union bound.
The reference ratio interval is [p1_lower/p0_upper,p1_upper/p0_lower], intersected
with the physically declared 0<kappa<=1. Zero denominators/lower endpoints and
empty intersections produce explicit conservative fallback, not clipping to a
positive ad hoc rate. Feed the ratio interval and retention intervals into the
previous sharp drift polytope. Save counts, intervals, endpoint witnesses and
primal/dual LP certificates. The four-interval image is conservative; no claim
of the shortest confidence interval or joint-set sharpness is made.

Primary gates at the largest reference budget: Wilson lower coverage bound >=.90
and median target-absence width <=.10 in every one of nine mechanism/ratio
cells. Report all smaller budgets without selecting an optimum after outcomes.
Reference event counts are not clock durations or monetary costs.

Defeat the measurement: repeat with q1/q0=.9 and assess (a) the incorrect common-
sensitivity assumption and (b) a correct q1/q0 in [.9,1] contract. Since the
observed reference ratio equals kappa*q1/q0, a bounded sensitivity ratio gives
kappa in [observed_lower/q_upper,observed_upper/q_lower]. Do not claim a false
measurement assumption becomes valid with larger sample size. Stress gates are
coverage only; report width and fallback frequency.

Before calling any biological cohort eligible, require initiated denominators,
retained counts in both intervention arms, a documented factor calibration,
independent loading-reference counts, and reference sensitivity/false-positive
information. Photobleaching-only or normalized survival data do not meet this
contract. Public GRID data DOI 10.5061/dryad.19st68k concerns CDX2 dissociation and
illumination intervals; its metadata does not supply this experiment. Its web
download returned HTTP 403 and official file API HTTP 401 in this environment;
no inaccessible outcomes were inspected, and no biological validation follows.
