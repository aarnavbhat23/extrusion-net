# Real trajectories do not separate the tested reversible subclass

The acquisition-date-held-out time-reversal screen **fails its detection gate**.
The primary unfiltered RAD21 baseline has adjusted p=0.224243; the filtered
sensitivity has p=0.763550. Neither supplies evidence for a new active mechanism.
Failure to reject is not evidence that chromatin is at equilibrium.

[Protocol](TRAJECTORY_REVERSIBILITY_PROTOCOL.md) committed at ac2578f;
[input amendment](TRAJECTORY_INPUT_AMENDMENT.md) at 882d0a2 before an evaluable
p-value. The initial run stopped on duplicate observation keys. Two unfiltered
tracks and one filtered track were excluded in full, without choosing a distance
among ambiguous localizations. One movie without an identifiable acquisition-date
prefix is excluded from both. All these exclusions remain in the result artifact.

| Baseline | Source rows | Development dates | Confirmation dates / movies / tracks | Adjusted p |
|---|---:|---:|---:|---:|
| unfiltered, primary | 28,040 | 7 | 14 / 23 / 138 | .2242431640625 |
| filtered, sensitivity | 17,605 | 7 | 14 / 23 / 139 | .7635498046875 |

The three lag features use exact frame separations of 20, 40 and 80 seconds.
Development tertiles in micrometres are [.303498653, .503330422] unfiltered and
[.287824961, .470224756667] filtered. Mean confirmation cycle scores are
[-.0113203, .0174551, .0174230] and [-.0168242, .0012831, -.00811984]. All
16,384 joint acquisition-date orientation flips are enumerated; the maximum
across lags gives the reported multiplicity-adjusted p-values. Tracks and frames
are not independent replicates. Standardization uses the conditional orientation-
randomization standard deviation sqrt(sum_date score_date squared); the frozen
zero-variance wording refers to this variance (all scores zero), not to empirical
sample variance across dates. Acquisition-date independence itself is an
assumption, not a verified biological-replicate design.

## What the conditional guarantee says

Let a latent stationary path Z have the same law as its reversal RZ. For an
emission kernel K that commutes with reversal,
P(Y in A) = integral K(A|z) dP(z) = P(RY in A).
Conditionally independent, time-homogeneous emissions satisfy this condition;
a reversal-equivariant missingness/selection rule also preserves it. A statistic
S with S(RY)=-S(Y) can therefore be tested by cluster orientation flips when
clusters are independent. Development thresholds must be independent of the
confirmation clusters; any fixed distance binning preserves the argument.

This is a consequence of established reversible-process theory. It does not
identify a motor state. An irreversible hidden cycle with a constant emission
has exactly the same observed path as a reversible constant process, so no
observable time-reversal diagnostic can universally identify hidden activity.
The tests include that counterexample. They also check a visible directed cycle,
an exact reversible chain, missing-frame reversal, full sign-orbit calibration,
and invariance to duplicating tracks inside a movie.

Observation rules matter: see [Hartich and Godec's comment](https://arxiv.org/abs/2112.08978)
and [the original authors' response](https://arxiv.org/abs/2202.02064) concerning
nonlocal time coarse graining. We use fixed-lag local binning and make no
entropy-production estimate. Tracking/QC and preparation can still violate the
combined null. The equilibrium subclass does not contain every BILD/switching
polymer model; no general separation from BILD was established or tested.

## Research decision

Stop this fixed three-bin diagnostic as a mechanism claim. Do not search new
bins/lags on these confirmation dates and relabel them independent confirmation.
A distinct learned biological rule still needs a matched predictive task against
active polymer and BILD-related controls, beyond a positive time asymmetry alone.
The next useful task is a perturbation-transport question: can a rule calibrated
on baseline and selected interventions predict a completely excluded intervention
at fixed genomic geometry, including localization noise and censoring? Freeze
that rule, its closest comparator, perturbation holdout and compute budget before
reading further condition outcomes. Require an independent study only after its
cadence, grouping and genomic differences are resolved. If no capability remains
after that comparator audit, retain the narrower identifiability work and pivot
rather than make a foundational architecture claim.

Source: Gabriele et al., DOI 10.5281/zenodo.5770531, CC BY 4.0. The unfiltered
file SHA256 is a67a887ac5ee7d13c637674486ad80a1efabab5e748a278761a0f482da181ee1.
The full source hashes, grouping, thresholds, track scores and test are in
results/trajectory_reversibility.json. The second independent trajectory study
has not been evaluated. This is neither independent biological replication nor
evidence that Extrusion Net surpasses existing learning architectures.

```bash
python study_trajectory_reversibility.py
python verify_trajectory_screen.py
python -m unittest discover -s tests -p 'test_trajectory_reversibility.py' -v
```
