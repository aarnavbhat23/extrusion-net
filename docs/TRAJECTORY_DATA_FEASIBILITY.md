# Public live-cell trajectory route: file access verified, mechanism not identified

The endpoint-order screen failed its power gate. A separate public-data search
now establishes that temporal observations can be obtained without first finding
a laboratory partner. This changes the data-access outlook; it does not turn
locus tracking into direct measurement of a cohesin motor's internal state.

Two independent releases were checked through their versioned Zenodo records,
then one small file per release was downloaded and verified against its supplied
MD5 and a recorded SHA256. The metadata-only reproducible inventory is in
results/trajectory_feasibility.json; raw files remain in ignored external/.
Both sources are CC BY 4.0 and retain the original authors' attribution.

| Source | Pilot file | Rows | Tracks | Grouping metadata |
|---|---|---:|---:|---|
| Gabriele et al., Dynamics of CTCF and cohesin mediated chromatin looping revealed by live-cell imaging | Rad21_0_hr.tagged_set.tsv | 17,605 | 183 | 34 movie_index values |
| Sabaté et al., Uniform dynamics of cohesin-mediated loop extrusion in living human cells | Adjacent_QC.csv | 12,269 | 150 | Track_pair; movie/batch mapping not yet resolved |

Gabriele release: https://doi.org/10.5281/zenodo.5770531 . The record documents
20-second frame spacing, 3D coordinates of both labeled loci, filtered/unfiltered
versions and several depletion conditions. This pilot downloads one untreated
RAD21 baseline, not all advertised conditions or a new validation cohort.

Sabaté release: https://doi.org/10.5281/zenodo.16949930 ; paper
https://doi.org/10.1038/s41588-025-02406-9 . Its pilot CSV includes both spots'
3D coordinates, distance, per-coordinate localization precision and a distance-
precision field. Frame indices do not by themselves establish cadence. Do not
infer biological-replicate IDs from the underscore structure of Track_pair
without checking the authors' mapping/code. An adjacent-locus pilot and a
mouse 515-kb baseline are not automatically comparable transfer cohorts.
The source paper already studies multi-locus kinetic consistency; claiming
that generic uniform kinetics is newly discovered would duplicate prior work.

No distance distribution, learned predictor, kinetic parameter, memory effect
or held-out outcome was evaluated in this metadata check. The other conditions
remain untouched by model selection. Completeness, censoring, tracking errors,
QC selection, cadence, movie/cell grouping and genomic span still need an audit.
The 2023 CTCF source workbook was also inspected at sheet/header level: it
contains encounter summaries, distributions and illustrative fluorescence traces,
not a verified complete ordered force/engagement event cohort.

The authors already developed Bayesian Inference of Looping Dynamics (BILD),
using a Rouse-polymer observation model to distinguish looping from proximity
([author-hosted paper](https://www.ashansenlab.com/uploads/1/2/5/3/125335035/ctcfcohesinloop_dynamics.pdf),
[implementation](https://pypi.org/project/bild/0.0.4/)). Inferring latent loop
states from these trajectories would therefore duplicate an established result.
A useful new candidate must survive that comparator and address a genuinely
untested predictive or intervention question.

## Next experiment worth freezing

Ask whether a restricted learned latent interaction rule predicts held-out
perturbation trajectory laws beyond passive polymer memory, tracking noise and
the authors' mechanistic/inference controls. Hold out whole movies/cells and
conditions, not frames; preserve an independent study for transfer only after
its acquisition and genomic differences are accounted for. Include temporal
subsampling so a supposed motor-memory signal must survive changes in cadence
rather than exploit measurement aliasing. The physical observation model must
include localization uncertainty and unresolved internal states.

This is a research question, not an announced architectural contribution or a
claim of identifiability. Predicting locus-distance history is insufficient to
identify a unique motor-state mechanism. Before fitting, specify which outcome
or perturbation could falsify the learned rule and the closest existing model.
A 9–10 novelty claim would require a distinct contribution surviving that
comparison and independent validation, not another name for state estimation.

```bash
python fetch_trajectory_pilot.py
```
