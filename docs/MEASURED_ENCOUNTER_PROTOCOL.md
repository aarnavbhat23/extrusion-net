# Measured encounter screen — frozen 2026-10-21

## Question and decision

Does an encounter rule conditioned on the current extrusion state predict
previously unseen experimental batches better than a force-only rule?
This is a bounded first screen of real single-molecule measurements, not a
new neural architecture, a causal intervention, or evidence for general
sequence learning. Existing force dependence is a published finding.

## Source and feasibility

Davidson, Barth, Zaczek et al., Nature (2023),
https://doi.org/10.1038/s41586-023-05961-5 . The authors' public analysis archive
https://doi.org/10.5281/zenodo.7409240 contains CTCF_overview.xlsx. Use its
encounter records; do not reconstruct raw trajectories from videos. Full raw
trajectories are available on request, not verified publicly here. Preserve
source attribution. Archive and workbook remain in ignored external/.

Metadata inspection before freezing found true Boolean Passed annotations,
force, elapsed frames, loop size, paths identifying batches, and repeated
molecule records. Some records aggregate repeated attempts. Do NOT expand
attempt counts into independent observations, use them as features, or infer
ordered history from worksheet order. Frames may use different experimental
cadences. State here means measured loop size and elapsed frames, not a
reconstructed latent history. Both variables may be measured retrospectively;
this is an observational prediction screen, not a prospective controller.

## Cohort and split

Primary: sheet CTCF, side N, exact Boolean Passed, finite nonnegative force,
elapsed frames and loop size; keep all such complete rows. No imputation.
Force must be >0 for the published exponential control. Group by the six-digit
experiment date in the path. The pooled 20220323_25 path is one group, never
split into purported days. Hold out every entire group once. Molecules do not
cross groups. Count repeat rows once each and give each molecule total weight
one, then each group equal total weight, in both fitting and primary scoring.
Report groups, molecules, rows, exclusion counts and force overlap. Secondary:
add side-C iCTCF rows and a side indicator using the same group split.

## Frozen controls and candidate

All fitted using training groups only; fixed hyperparameters, no outcome-driven
selection. Pass probability:
1. smoothed constant;
2. author's functional family exp(-F/F0), F0 fitted by likelihood;
3. logistic with log(1+F/.025);
4. logistic cubic in log(1+F/.025), ridge penalty 1 per effective row.
Secondary controls add a C-side offset (including exponential via an F0 ratio).

Candidate: same linear-force logistic plus log(1+elapsed/100), loop size/10,
and force-by-loop-size interaction. Ridge penalty 1, training-only feature
standardization, unpenalized intercept. The force-only cubic guards against
mistaking a flexible force curve for state dependence. Report a same-input
histogram gradient boosting control (max_iter=100,max_leaf_nodes=7,
learning_rate=.05,min_samples_leaf=10,l2_regularization=1,random_state=21);
this tests whether a claimed benefit exceeds a generic predictor. This generic
control is included in the primary best-control comparison.

Primary metric: equal-group, molecule-weighted held-out log loss. Secondary:
Brier score and individual group changes. Bootstrap whole held-out groups
10,000 times, seed 21, paired candidate minus each control. Success requires
>=10% reduction relative to the best mean control, lower log loss in a strict
majority of groups, and a 95% paired group-bootstrap interval for improvement
above zero against every control. The sample is small; this screen cannot
establish population or causal significance. Do not modify the gate afterward.
No gate passage means stop expansion of this candidate, retain an honest
measured-data benchmark, and require ordered trajectories before pursuing
history learning. No unseen test data remain after this screen.

## Novelty boundary

Force-sensitive barriers are Davidson et al.'s result. Dynamic barriers and
lifetime effects already appear in Rahmaninejad et al. (2025),
https://doi.org/10.1101/gr.280108.124 . Logistic interaction models and grouped
validation are standard tools. A successful screen would motivate a future
learned kinetic rule; it would not itself establish a novel modeling paradigm.
