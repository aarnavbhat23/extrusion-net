# Joint testing removes the old witnesses, but not finite-data ambiguity

**All twelve rectangular-calibration witnesses are excluded by the new joint
likelihood set. New feasible pairs nevertheless differ by 2.12–5.64 percentage
points in intervention absence across all six datasets.** Every pair exceeds
the original two-point precision requirement. This stops the claim that swapping
calibration statistics alone has solved prediction in this experiment.

The [protocol](JOINT_LIKELIHOOD_PROTOCOL.md) was committed at `36acd46`, before
evaluating this reanalysis. These are the same six existing datasets. A native
boolean serialization repair is recorded in [the run note](JOINT_LIKELIHOOD_RUN_NOTE.md);
the identical rerun reproduces every printed outcome. There was no split-seed
selection or replacement of failed cases.

| mu / pL / pR / h | Repeat | Old witnesses retained | New attained span (points) | Generating log-likelihood slack |
|---|---:|---:|---:|---:|
| .21 / .3 / .7 / .4 | 0 | 0/2 | 4.320 | 4.004 |
| .21 / .3 / .7 / .4 | 1 | 0/2 | 5.637 | 7.781 |
| .45 / .65 / .35 / .8 | 0 | 0/2 | 3.234 | 1.027 |
| .45 / .65 / .35 / .8 | 1 | 0/2 | 4.548 | 10.031 |
| .8 / .4 / .6 / .55 | 0 | 0/2 | 2.121 | 3.009 |
| .8 / .4 / .6 / .55 | 1 | 0/2 | 4.302 | 9.962 |

## What the guarantee covers

Each original 10,000-episode assay is split into independent 5,000-episode
training and evaluation samples by a seeded random partition of its counts.
Training jointly fits both gamma=.3/.1 assays, without using generating
parameters. The held-out likelihood ratio defines the established split
likelihood confidence set at alpha=.05. Conditional on the model and independent
sampling, a prospectively fixed procedure covers the parameter with probability
at least .95 per dataset. This is
[universal inference](https://arxiv.org/abs/1912.11436), not a new statistical method.
The post-outcome choice to reanalyse these data is not prospective confirmation.
The mathematical statement averages over the random partition and the data;
it is not a coverage claim conditional on the recorded aggregate counts or a
particular pseudorandom seed.

All 18 training starts converged. 47/48 local endpoint searches report success;
all statuses are saved. Acceptance depends on directly checked feasibility, not
optimizer status. The reported spans are lower bounds on set diameter, never
confidence interval widths or global extrema. All six generating parameters
are retained, which is a replay check, not validation of nominal coverage.

The old and new sets are not nested. The old bounds use all counts, while the
new test spends half on fitting and avoids marginal Bonferroni rectangles.
Therefore the width comparison is not a pure measure of conservatism, and no
universal information lower bound follows. Neither arbitrary misspecification
nor biological observation noise is covered.

## Structural versus practical identification

A new [explicit inverse](CAPTURE_IDENTIFIABILITY.md) proves that a single exact
snapshot distribution already identifies this family, except for an irrelevant
parameter when capture is absent. Its informative beyond-barrier states have
only **1.95–21.08 expected observations per 10,000-episode assay** in these cases;
the observed range is **0–23**. This explains a potential conditioning problem,
without proving these statistics are the most efficient possible estimator.

The inverse passes 22 independent 70-digit forward/inverse checks, including
one-sided capture, no capture, complete capture and changed removal ordering.
The proof assumes known loading, stepping, barriers and observation clock. It
does not identify a biological mechanism from ordinary contact maps.

## Reproduction and decision

Run `python study_joint_likelihood.py`, `python verify_joint_likelihood.py`, and
`python study_capture_identifiability.py`. The saved
[likelihood artifact](../results/joint_likelihood.json) contains every split,
fit and witness. A separate decimal forward implementation verifies all saved
likelihood slacks, intervention spans, input hashes and split reconstruction at
70 digits. The [inverse audit](../results/capture_identifiability.json) records
every boundary check and informative-state count.

Do not add further hidden states to address a precision failure in an already
structurally identifiable model. Next distinguish measurement efficiency from
family misspecification: compare accessible observation schedules and a strong
full-data inference control, then require prediction of an excluded measured
perturbation against the closest published kinetic model. A generic latent-state
model is already represented in the literature, including
[Shah and Tortora et al.](https://doi.org/10.1016/j.molcel.2026.03.037).
A useful new mechanism or a scalable, sharper guarantee remains to be established.
