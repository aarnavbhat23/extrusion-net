# Frozen screen: a falsifiable equilibrium subclass

Scope: this is a diagnostic prerequisite for a mechanistic model, not a new
architecture and not a comparison against all of BILD. Freeze before accessing
outcome distributions in the new unfiltered file. The previously downloaded
filtered RAD21 baseline has only undergone a metadata feasibility check.

Source: Gabriele et al., DOI 10.5281/zenodo.5770531, CC BY 4.0.
Primary file Rad21_0_hr.unfiltered.tagged_set.tsv, source MD5
7e9a1063b1035b5c380a0a45c067a311; filtered baseline is sensitivity only.
Twenty seconds per frame and micrometre distances follow the source README.

Null: each acquisition-date cluster has a path distribution invariant under
joint reversal of all its tracks (including observation/missingness patterns).
Different dates are assumed independent, not different frames or tracks within
a movie. Group movies by a leading YYYY_MM_DD identifier; exclude unidentified
movie dates with an explicit count. Dates are acquisition clusters, not verified
independent biological replicates. No causal interpretation follows.

Order dates by SHA256('reversibility-v1:' + date); the first ceil(D/3) dates are
development, the remainder confirmation. Fit only two pooled development-distance
tertiles. Apply these unchanged to confirmation tracks. Use exact frame-lag pairs
at 1, 2 and 4 frames, without interpolating missing frames. For each ordered pair
of distinct bins, score +1 on 0->1, 1->2, 2->0, -1 on their reverses, 0 otherwise.
Average within track, within movie, then within acquisition date with equal
weights at each level. Reverse all tracks in a date together.

Primary test is the maximum absolute standardized mean across the three lags.
Enumerate every date sign flip if D_confirm <= 20; otherwise use 65,535 seeded
sign flips (seed 2117), with the Monte Carlo plus-one correction. Zero-variance
features score zero. Familywise alpha=.05 over these three lags. Require at least
8 confirmation dates; otherwise report insufficient clusters, not a pass.
Sensitivity: run the same pipeline separately on the filtered baseline, without
using its outcomes to change the primary test. No additional conditions or new
features will be selected following results.

Separation guarantee: a stationary reversible latent process observed through
a time-homogeneous conditionally independent emission kernel remains reversible.
A time-reversal-equivariant selection rule preserves the property. Thus rejection
can falsify that *combined* model/measurement subclass. It cannot distinguish
active extrusion from nonstationary preparation, tracking directionality,
non-equivariant QC or another nonequilibrium process. Failure to reject cannot
establish equilibrium: an irreversible hidden cycle can be invisible in distance.
This uses established time-reversal/detailed-balance principles, not a new theorem
of statistical mechanics. Full BILD and switching polymer models are not all in
this null; there is no general BILD separation claim.

Positive and negative controls: independent stationary three-state reversible
chains and directed cycles, including a hidden directed cycle with a constant
emission. Unit tests must check exact score reversal, missing-frame handling,
cluster duplication invariance, exact randomization calibration on a small orbit,
and blindness to the constant-emission alternative. These are implementation
checks, not independent biological confirmation.

Success for this screen: primary adjusted p <= .05 with sufficient clusters.
Even a pass only earns investigation of mechanism and acquisition artifacts.
A 9-level contribution would additionally need a distinct update law, separation
from the closest active/polymer models, predictive perturbation validation, and
an independent, appropriately matched study.

Prior art: Martinez et al., Inferring broken detailed balance in the absence of
observable currents, Nature Communications 10, 3542 (2019),
https://doi.org/10.1038/s41467-019-11051-w ; Gabriele et al.,
https://doi.org/10.1126/science.abn6583 .
