# Measured iFRAP transfer screen — frozen 2026-10-29

Freeze before fitting or scoring the curated source curves. Source CSV SHA-256:
`5f9e49088e2218d6f35bca61f3039da72864906de711740437949e11cb47ec0b`.
Use the dataset-release-v0 day partition unchanged: 15 development, six reserved
test and eight quarantined movies. No quarantine relabeling, no nucleus/frame
splits and no test-driven hyperparameter selection. All five conditions with
eligible development/test days are reported. Other conditions remain development
only. Test days are technical replication, not independent biological cohorts.

## Observable and candidate mechanism

For each movie, average source unbleached/bleached ROI intensities across nuclei.
Let D(t) be their difference, and C(t) the shared unbleached-control minus
background. Relative to source frame 5, define X=D/D5, Z=C/C5 and Y=X/Z.
Time is in minutes since frame 5. The existing normalized source curve equals Y.
Use all 245 postnormalization frames without clipping negative values.

A restricted reaction-limited measurement model is X(t)=Z(t)S(t)+noise, where
S is a nonnegative exponential mixture. This does NOT identify loop extrusion,
biological removal, diffusion, binding states, or the physical origin of its
rate components. Measured controls are conditioned upon, not known noise-free
biology. No claim of unconditional future-image forecasting is allowed.

Fit equal-weight development movies within each condition, retaining all starts:
1. one exponential exp(-kt);
2. plateau plus exponential p+(1-p)exp(-kt);
3. two exponentials w exp(-k_s t)+(1-w)exp(-(k_s+gap)t);
4. the same two-exponential family fitted to X-Z*S (raw-contrast objective);
5. the empirical condition-mean development curve (strong nonparametric control).

Rate bounds [1e-6,5] per minute; weight/plateau [0,1]. Fixed starting rates
[.001,.01,.1] for single/plateau and slow rates [.001,.01,.1], gaps [.01,.1],
weights [.25,.75] for mixtures. Least-squares, at most 1500 evaluations/start.
Retain optimizer failures. Select the lowest objective among converged starts;
if none converge, report failure rather than an omitted movie.

Implementation clarification recorded 2026-10-30: the mixture bounds apply to
the optimized slow-rate and gap coordinates separately; their sum can in
principle reach 10/minute. All selected nominal mixture fast rates are below
2.22/minute. No outcomes were used to change these implemented bounds.

The raw-contrast objective equals Z^2-weighted squared error on Y. It is a
measurement-aware control, not a new architecture or inference principle.
Verify this identity numerically. Its interpretation assumes the numerator
noise model; noisy denominator/reference values can still bias estimates.

## Frozen comparisons and gates

Primary outcome: normalized-curve MSE, reported per held-out movie and averaged
equally across the five conditions. Also report MAE and the worst movie.
No frame-level p values or independence claims. Do not select a winning baseline
per test movie: compare the candidate separately against every declared control.
Usefulness screen: raw-contrast fitting reduces condition-balanced MSE by at
least 10% against BOTH the ordinary mixture and empirical mean, and improves
at least four of five conditions against each. This is a computational screen,
not a biologically approved effect size or statistical confirmation.

Development-only sensitivity: impose a +/-5% linear drift on Z over the recorded
window, leave X unchanged, refit the two mixture objectives, and record rate/curve
changes. This is a deliberately imposed error, not an estimated instrument-error
distribution. There are only five held-out conditions; no paper-level superiority
or confidence claim follows from passing this small screen.

Independent checks: original-source reaggregation, saved-parameter prediction
replay, split identity, exact weighted-loss equivalence and synthetic recovery.
If controls match or beat the candidate, stop this as a novel-method claim.
Keep the earlier internal-calibration theory separate: this assay does not
measure its required state-specific hazards or calibrated perturbation strength.
