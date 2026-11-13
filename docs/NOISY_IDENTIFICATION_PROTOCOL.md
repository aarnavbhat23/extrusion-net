# Noisy contact-map stress protocol

Written before this study's results. This is an exploratory synthetic experiment,
not an externally registered confirmation or biological validation.

Question: does intervention selection along a fitted observational equivalence
family still improve prediction when the baseline contact map is noisy?

Use 10 independent mechanisms at each of 8, 16 and 32 sites (matrix dimensions).
Loads are uniform [0.4,0.9], directional gates uniform [0.4,0.7]. Add independent
Gaussian noise to every observed upper-triangular contact, including the diagonal,
with absolute standard deviations 0.001, 0.01 and 0.05. Intervention observations
have the same noise standard deviation. Seed 20261026; no model tuning on results.

Compare raw three-band algebra with constrained least-squares projection of all
contact bands. Fit loads in [0.05,2] and gates in [0.005,0.995] from two fixed starts;
retain the lower objective. The bounds are numerical modeling assumptions, not
known population bounds. Record convergence, iterations and baseline residuals.

Select a single directional gate to cut using the fitted family's response span,
a fixed boundary, or a uniform random site. All strategies share the same scalar
measurement noise. Fit the baseline map jointly with the chosen intervention.
Also report the midpoint of the projected equivalence family without intervention.
Evaluate other local intervention responses, excluding the selected site. Report
paired differences by mechanism, including failures, not only favorable seeds.

All computations are synthetic. Selection sees the noisy baseline only. Truth is
used only for generating observations and scoring. Do not claim confidence bounds
from a projected point estimate, a new architecture, or novelty priority.

## Exploratory extension after the first stress test

The span rule was inconsistent under noisy baselines. Add a local linear design
criterion using the full baseline Jacobian and its one-dimensional null space.
For candidate intervention gradient g, null direction n and P=(J'J)^+, the
predicted covariance of response gradient h is
sigma^2 [hPh' - 2a hPg' + a^2(1+gPg')], a=(hn)/(gn).
Choose the site minimizing the mean predicted variance of all other local
interventions. This is classical local optimal experimental design specialized to
the routing model, not a novelty claim. Use 10 fresh mechanisms per size with seed
20261028, the same noise grid, all four selection rules, and both held-out and
all-site errors. Retain the first study separately. The linear approximation may
fail near constraints or at high noise; report those failures.

## Robust design extension after the local-design experiment

The single-point local criterion also failed in some conditions, notably size 16
at noise 0.01. It is evaluated at an arbitrary member of an observationally
indistinguishable family. Test a worst-case local variance criterion over nine
fixed interior fractions (0.05 through 0.95) of the admissible family. This is a
grid approximation, not a certified continuous minimax solution. Near-singular
families and nonlinear response errors remain possible.

Use fresh seed 20261030, the same 10 mechanisms per size/noise grid and all five
rules. For every rule, jointly fit from three equivalence-family starts, retaining
the lowest observed residual objective, to check optimizer sensitivity fairly.
This is another exploratory extension, not an independent replication of a
preselected winner. Do not select a favorable seed or suppress earlier results.
