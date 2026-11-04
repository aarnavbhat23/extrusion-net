# Measured iFRAP transfer: useful audit, failed method gate

The frozen [protocol](MEASURED_IFRAP_PROTOCOL.md) was committed as `ddf9a06`
before fitting. Source curation uses 15 development movies, six reserved movies
across four acquisition days and five eligible conditions, and eight quarantined
movies. All recorded outcomes are retained. Test frames and nuclei are not
treated as independent observations. Source labels/date ambiguities are unresolved.

| Method | Condition-balanced held-out MSE | Worst-movie MSE |
|---|---:|---:|
| Empirical development mean | .023291 | .062435 |
| Single exponential | .024583 | .066863 |
| Plateau plus exponential | **.019975** | .061107 |
| Two-exponential mixture | .020272 | .061999 |
| Raw-contrast mixture | .020047 | **.055717** |

The raw-contrast candidate reduces MSE **13.93%** against the empirical mean,
improving all five conditions. Against the ordinary mixture, its reduction is
only **1.11%**, improving **two of five** conditions. It does not beat the plateau
control on aggregate MSE. The prespecified 10% / four-of-five gate against BOTH
mixture and empirical controls therefore **fails**. The slightly better worst
movie is a secondary diagnostic, not a replacement success criterion.

There are 210 nominal optimizer starts and 336 development-only sensitivity
starts; all converge. Imposed +/-5% reference drift changes fitted development
curves by up to .03428 in normalized signal. Those perturbations are synthetic,
not measured instrument error. The raw-contrast loss equals conventional
Z-squared-weighted normalized squared error; maximum numerical discrepancy is
6.94e-17. There is no new architecture or statistically demonstrated superiority.

## What the measured data establish

We can reproduce source fluorescence processing and predict normalized curves
on reserved technical days. This checks a real observable rather than transferring
an unexamined simulator label onto biology. It does not validate loop extrusion
states, biological removal rates, calibrated intervention effects or a universal
learning mechanism. A reaction-limited exponential family remains an assumption;
reference dynamics and diffusion require independent checks.

The [measurement-operator deduction](MEASURED_OBSERVATION_THEORY.md) gives an
explicit obstruction: two nonnegative optical/biological/reference rate triples
produce identical complete raw fluorescence trajectories but different responses
to a selective intervention. An exact-data target interval has width about
8.85 percentage points in the stated toy example. This is a restricted
structural result and source of experimental requirements, not empirical evidence
that the source assay has those particular rates or a priority claim.

## Verification and reproduction

- [Complete fits and sensitivity starts](../results/measured_ifrap.json).
- [Per-movie outcomes](../results/measured_ifrap_movies.csv) and
  [all 7,350 saved prediction rows](../results/measured_ifrap_predictions.csv).
- [Independent source reaggregation and prediction replay](../results/measured_ifrap_verification.json): maximum discrepancy 1.28e-15, without refitting.
- Four new tests cover known-case recovery, weighted-loss equivalence, model
  limits and exact raw-signal/intervention ambiguity. The full local suite has
  114 passing tests.

The [source-only acquisition package](../datasets/ifr_v01/README.md) downloads
the original public CSV and reconstructs the curated split:

```sh
python datasets/ifr_v01/build_extrusion.py
python study_measured_ifrap.py --dataset-root datasets/ifr_v01
python verify_measured_ifrap.py --dataset-root datasets/ifr_v01
```

## Paper decision

Defer submission and stop this weighted-fitting variant as a novel-method claim.
The appropriate direction remains a focused identifiability/experimental-design
study: establish useful finite-sample predictions when perturbation strength and
observation loss are uncertain, with a measured observation operator. A repeat
measurement that independently constrains reference dynamics is more valuable
than additional curve-model variants. The earlier internal-calibration result
requires state-specific rates that this dataset does not provide.

Research quality improves through a real-data test, strong controls, preserved
failures and an explicit measurement limitation. Originality has not earned a
higher rating. Independent biological validation and closest-prior-work separation
remain unmet submission gates.
