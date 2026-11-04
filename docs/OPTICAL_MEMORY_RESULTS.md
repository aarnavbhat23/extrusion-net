# Optical attribution: a sharper boundary, another failed order gate

Protocol `f006846` and estimator/tests `dc22e8c` were committed before the larger
screen. This implements the six-step counterexample → measurement → adversary →
identification → frozen comparison → prior-work workflow. All counts are
synthetic; no new biological observation is claimed.

## What happened

The two-state biological and optical models have exactly equal primary
observation laws. All 340 enumerated short pulse/observation strings agree,
and an independent hidden-path expansion checks the forward calculation.
The [proof](OPTICAL_MEMORY_THEORY.md) covers arbitrary and adaptive sequences.
Primary order effects therefore do not identify biological memory.

With an independently selective calibrated reporter, the order contrast is
identifiable under the stated calibration/error interval. But the frozen power
gate requires at least 80% power in four of six cases at n=2000 per arm.
It reaches that threshold in **three of six**, so the gate **fails**.

| Capture a | Release b | Reporter order power | Direct A-versus-baseline power |
|---|---:|---:|---:|
| .2 | .3 | 0 | 1 |
| .2 | 1.2 | .6915 | 1 |
| .6 | .3 | .4725 | 1 |
| .6 | 1.2 | 1 | 1 |
| 1.5 | .3 | .999 | 1 |
| 1.5 | 1.2 | 1 | 1 |

The simpler control tests presence of binding, not the same order-effect target.
It shows that pulse ordering is unnecessary for separating these particular
models once a trusted state reporter is available. Equal endpoint counts do not
establish equal reporter/instrument costs. Rates are dimensionless synthetic
choices, not estimated biological constants.

Across all in-contract reporter settings, minimum observed interval coverage
is 99.9%; maximum optical-null rejection is .05% (Wilson upper .283%). These
are conservative fixed-time intervals with assumed calibration bounds, not
calibrated real-assay coverage. Each rate uses 2,000 experiment replicates.

When the optical state is allowed to affect the added reporter in the same way
as binding, the guarantee no longer applies. In three n=2000 settings the method
attributes binding in all simulated optical-only experiments. That is an
explicit out-of-contract failure, not a contradiction of the conditional bound.
The full-law counterexample explains why no estimator can solve that expanded
model class without a new restriction or observation.

## Reproduction and originality

```sh
python study_optical_memory.py
python verify_optical_memory.py
python -m unittest discover -s tests -v
```

- [Complete outcomes](../results/optical_memory.json), [summary](../results/optical_memory_summary.csv), and compressed original count arrays are retained.
- [Independent replay](../results/optical_memory_verification.json) checks all 168,000 experiment replicates without importing the study estimator.
- 118 local tests pass, including four new independent solution/bound tests.
- [Equation/protocol prior-work audit](OPTICAL_MEMORY_PRIOR_ART.md) identifies substantial overlap with HMM equivalence, binding/photophysics inference and experimental design.

## Research decision

Stop pulse-order novelty/efficiency claims for this prototype. Retain the
counterexample and calibration requirements. Next choose a physically achievable
state or optical reference, quantify its cross-talk and perturbation specificity,
and derive which intervention predictions it identifies. A generic extra channel
or another pulse sequence is not sufficient. Novelty stays approximately 6.5/10;
research quality approximately 7.5/10, with roughly one-point judgment uncertainty.
No manuscript submission or new architecture claim is justified.
