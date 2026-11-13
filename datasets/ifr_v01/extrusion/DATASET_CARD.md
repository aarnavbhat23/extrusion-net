# Extrusion observation dataset — development release 0.1

## What this is

A curation of measured RAD21 iFRAP fluorescence observations from the public
[Louafi and Giorgetti dataset, version 1](https://doi.org/10.17632/spp4x8jpdy.1).
The experiments concern cohesin cofactors. The source table contains 30,750
measurements from 123 nuclei in 29 movies and seven source-defined conditions.
These are **not** 30,750 independent samples, individual molecular residence
times, loop lengths, extrusion trajectories, or ground-truth removal events.

Our additions are persistent observation identities, trace/movie/day grouping,
explicit time units, source-quality flags, reproducible source normalization,
safe default loaders, and an evaluation partition. We did not collect the
source experiments or create new biological labels.

## Source and processing

- Original file: `combined_dataframe.csv`; publisher SHA-256 is recorded in
  `source_manifest.json` and enforced before every build.
- Source data license: CC BY 4.0. Credit Nessim Louafi and Luca Giorgetti; retain
  the DOI and indicate that this is a derived curation.
- Processing reference: `NessLfy/cohesin_residence_time`, commit
  `a2897e6da9b7074c88770343348f7ed308e24043`. The BSD 3-Clause notice is retained
  in `../source/ifr_LICENSE` for accompanying source material.
- The instrument description specifies 250 frames, one every 30 seconds, and
  bleaching after four frames. We expose time from the first frame, **not** an
  invented exact bleaching timestamp. The notebook's plot uses minutes; one
  inline comment incorrectly calls a half-minute increment seconds.
- The word `interp` in the source columns concerns interpolation of ROI
  positions between manual updates. Intensities are then sampled from images
  at each frame. The table is processed measured fluorescence, not a set of
  fabricated intermediate intensity measurements. Its rows remain correlated.
- The source pipeline selected/excluded nuclei before this table was produced.
  This release cannot recover discarded nuclei, an exclusion denominator, or
  the causes of signal disappearance.

## Tables

`observations.parquet` preserves all nine original columns and their values.
`source_row_id` is the zero-based original data-row index. `movie_id` is a stable
hash of the complete source filename. `trace_id` combines movie and nucleus.
The source `nucleus` identifier alone must never be used as a global key.
`time_since_first_frame_s = 30 * frame`. Fluorescence is in source instrument
intensity units, not molecule counts.

`movie_curves.parquet` has 250 rows per movie and replays the authors' operation:
average unbleached and bleached intensities across that movie's nuclei; subtract;
divide by the shared unbleached-control-minus-background signal normalized to
its first four frames; then normalize to frame 5. All frames are retained, and
`source_postnormalization_window` identifies the authors' frame-5 onward window.
`time_since_normalization_s` is relative to frame 5, not a biological event.
Frame-5 normalization must not be used to claim forecasting before frame 5.

`movies.csv`, `nuclei.csv`, `source_issues.csv`, `condition_coverage.csv` and
`splits.json` make provenance, grouping, ambiguity and coverage inspectable.
No biological rate estimate is embedded as ground truth.

## Source-quality issues and partition

Four movies have treatment labels inconsistent with filename hints: both WAPL
movies dated 2024-04-26 and both PDS5 movies dated 2024-05-27. We cannot determine
whether the filename or table label is correct. Both are preserved without
relabeling, and those days are quarantined from evaluation.

Two filenames start with malformed `202402202`. Their dates are left missing.
The valid `20240220` group is also quarantined because it could share the same
underlying day as that malformed identifier. This is conservative grouping,
not a declaration that the valid filenames are wrong.

| Partition | Movies | Nuclei | Observations |
|---|---:|---:|---:|
| Development | 15 | 65 | 16,250 |
| Reserved test | 6 | 21 | 5,250 |
| Quarantine | 8 | 37 | 9,250 |

Whole acquisition days remain together. Partitioning used metadata and condition
coverage before predictive fitting; it is neither random by row nor chronological.
No downstream model has been evaluated during this dataset build. The public
reserved partition is not a securely blinded test. No independent donor or
laboratory replicate is claimed. Days inferred from filenames are technical
grouping identifiers, not proof of biological independence.

Five conditions appear in development and reserved test. `WAPL | - aux` and
`PDS5 | - aux` lack a clean second acquisition day, so **they have no valid
day-level holdout here**. There is no universal calibration partition. Methods
needing calibration must declare an eligible grouping strategy before use;
do not tune on reserved test days or split nuclei to manufacture replication.

## Valid next experiments

1. Specify how the biological model generates the measured iFRAP observable,
   including bleaching and the ROI/control normalization. Population recovery
   cannot automatically be treated as individual-molecule survival.
2. On development days, compare a conventional fluorescence-kinetic model with
   observation-aware inference under a shared likelihood and observation operator.
   Retain uncertainty at the movie/day level.
3. Freeze estimands, fitting procedure and failure criteria; evaluate only the
   five conditions with eligible reserved-day coverage. A day-held-out result
   would support technical transfer, not causal identification of removal.
4. Resolve the label/date discrepancies using authoritative source records or
   an explicitly authorized author query. Do not silently substitute filenames.
5. To validate removal-versus-observation separation, obtain repeated observation
   schedules, calibrated perturbation strengths and an independent physical
   reference. Those labels are **absent** from this release.

This dataset improves reproducibility and exposes falsifiable assay questions.
It does not establish a new learning architecture or make a paper ready by itself.
