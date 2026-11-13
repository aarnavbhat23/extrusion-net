# Counted-reference audit and measurement contract — 2026-11-03

Scope: this note audits the 9/286 counterexample in [COUNTED_LOADING_THEORY.md](COUNTED_LOADING_THEORY.md)
and the loading-drift sharpness proof in [LOADING_TARGET_THEORY.md](LOADING_TARGET_THEORY.md) as of
`d9cd4d3`. It builds a measurement-contract table for one concrete assay and chooses one next
experiment. The usefulness rule (§2) and the pilot go/no-go criteria (§4) are frozen by this commit,
before any pilot data exist. No new biological data were obtained. Everything below is exact-data
algebra or a design proposal. None of it is validation.

## Bottom line

1. **The counterexample and the sharpness proof are correct within their stated scope.** Independent
   exact-rational checks found no error. Three statements of scope need correcting; they are listed below.
2. **The current contract cannot give a useful prediction, however many counts are collected.**
   With f known to ±10% and q1/q0 in [0.9,1], the exact-data width (infinite counts) is already
   **6.26 pp** for mechanism M1. The true target is only **2.60 pp** from a naive linear-in-dose
   extrapolation and **4.55 pp** from "no extra effect beyond f". To get below 2.6 pp, the factor f,
   the target multiplier g and the arm loading ratio must all be calibrated to roughly **±1–2%**.
3. **In the assay most likely to meet the contract, the counted reference is not needed.** In a
   single-molecule reconstitution where the release factor is added after a loading pulse, κ=1 and
   q1/q0=1 hold by design. A zero-dose arm turns f into a linearity check. The target then reduces to
   a classical Stern–Volmer extrapolation in dose. Direct dwell-time imaging is a stronger control:
   it separates the two counterexample worlds without any reference.
4. **The next experiment should be an assay-feasibility pilot with go/no-go criteria frozen in
   advance, not another simulation.** The open gaps are biological: equal susceptibility under the
   release factor, dose linearity at the 1–2% level, and which event counts as "loss". Simulation
   cannot settle these. Running the pilot needs a lab collaborator. Until it passes, keep the paper
   to a narrow theory/design note.

## 1. Audit

### Replay of the existing record
In a scratch copy of the checkout, `verify_counted_loading.py` regenerates 7,200 reference
datasets, 10,800 interval images and 21,600 certificates. Its output is byte-identical to the
committed verification JSON. `verify_loading_research.py` also passes, and so do the 13 unit tests
in `test_counted_loading.py` and `test_loading_target.py`. The full 135-test suite was **not**
rerun.

### Independent exact checks (`audit_counted_contract.py`, no repo imports, no optimizer)
| Check | Result |
|---|---|
| Counterexample: both worlds physical, four probabilities equal (9/28, 1/5, 27/40, 27/50), targets 35/44 vs 43/52 | Exact; gap 9/286 = 3.147 pp |
| Are A and B special? | They are exactly the two endpoints of the exact identified set when f=2 is known and q1/q0∈[0.9,1]. So 9/286 is the full structural width of that contract. |
| Necessity: 4,000 random physical mechanisms map into the drift polytope, with target = 1/(U+W) | All pass |
| Sufficiency: 4,000 feasible points (394 on the W=0 face, others on the J and W bounds) reconstruct to physical mechanisms that reproduce R0, R1 and the target exactly | All pass |
| The closed-form exact-data interval equals exact LP vertex enumeration | 40/40 |

### Corrections to the stated scope
- **κ≤1 is stronger than the proof needs.** The construction sets ℓ=ℓmax and uses only
  `kU·ℓmax ≤ 1`. If that condition fails, the LP is still a valid outer image, but the sharp set
  gains the bilinear constraint `U·V ≥ J`. The audit includes an exact example with κ=5/4 and ℓmax=4/5.
- **κ≤1 does identifying work near κ=1.** When true κ=1, the [0.9,1] sensitivity contract adds zero
  width only because κ≤1 is imposed. Allowing κ up to 1/ℓmax (ℓmax=0.95) gives 0.85–1.62 pp instead.
  If a treatment can raise loading even slightly, that assumption must be stated and justified.
- **The target multiplier g is treated as exact while f is uncertain.** For a measured, dose-defined
  target, g carries the same calibration error as f. With equal relative error, the f-induced width
  grows about 1.7× (M1, ±10%: 3.02 → 5.23 pp). The LP still extends linearly in M=μ/(ℓγ); the audit
  uses that closed form.
- **"Complete observation-law equivalence" holds for presence-only readouts.** If the retained motor
  state is also recorded and its generator is unknown and reversible, the counterexample still
  survives (exact two-state matching world). If the state dynamics are irreversible or independently
  calibrated, and the perturbation leaves them unchanged, normalized state shapes identify c1/c0. κ
  then follows from R1/R0 with no reference (the audit recovers κ=4/5 exactly).
- **False positives in the reference are second order at small rates.** At a 5% false-positive rate
  among unloaded episodes and equal sensitivity, the same four probabilities are reproduced with
  κ=98/125 instead of 4/5. The target moves by 0.69 pp. The naive ratio is biased, but not by much
  next to the current widths.
- **Reference exchangeability is an extra assumption.** If reference and retention cohorts load
  differently within an arm, the ratio becomes κ_ref·q1/q0. This enters exactly like a
  sensitivity-ratio error and belongs in the contract.

## 2. Structural floor and calibration requirements

Exact-data widths in pp (infinite counts), for the frozen counted-study mechanisms with κ=0.8.
Usefulness reference: true target vs naive extrapolations, computed at κ=1.

| Mechanism (μ, δ, ℓ) | f ±10%, κ exact | q1/q0∈[0.9,1], f exact | both (current contract) | gap to "no extra effect" | gap to linear-in-dose |
|---|---:|---:|---:|---:|---:|
| M1 (1/5, 1/5, 3/4) | 3.02 | 3.15 | **6.26** | 4.55 | **2.60** |
| M2 (2/5, 1/10, 17/20) | 3.22 | 2.28 | 3.89 | 5.31 | 5.31 |
| M3 (13/20, 3/25, 19/20) | 2.44 | 1.23 | 2.44 | 4.54 | 5.52 |

What calibration buys in M1, the binding case (true κ=1, two-sided arm-loading uncertainty, g with
the same relative error as f):

| f and g tolerance | arm loading ±1% | ±2% | ±3% | ±5% |
|---|---:|---:|---:|---:|
| ±1% | 1.19 | 1.86 | — | — |
| ±2% | 1.71 | 2.38 | 3.05 | 4.39 |
| ±4% | — | — | 4.08 | 5.41 |

These floors come before any sampling error. A prediction that can reject linear extrapolation in M1
therefore needs f, g and the arm loading ratio all known to about ±1–2%. That is an explicit,
demanding assay requirement. The present ±10% contracts cannot meet it.

**Usefulness criterion (frozen).** The 95% target interval must exclude
both naive extrapolations, both computed from the same calibration data: (a) target absence equals the
dose-2 absence; (b) absence is linear in dose. This is operational, needs no target data, and is not
tuned on outcomes. For the frozen mechanisms it implies widths below 2.60 / 5.31 / 4.54 pp.

## 3. Measurement-contract table: one candidate assay

**Assay:** single-molecule reconstitution on surface-tethered DNA in a TIRF flow cell. Purified,
fluorescently labeled human cohesin is loaded with NIPBL in a short pulse, then washed. The loss knob
(WAPL–PDS5) is then added at dose d∈{0,1,2}, with d=3 held out as the target. Each field is imaged
once after a dark interval drawn per chamber. Reference chambers are imaged immediately after the
wash. This is the assay class of the human-cohesin loop-extrusion experiments in [PRIOR_WORK.md](PRIOR_WORK.md).
Wutz et al. 2026 ([DOI](https://doi.org/10.1016/j.molcel.2026.03.025), found via PubMed) report from
single-molecule imaging that PDS5 stops loop extrusion by promoting NIPBL dissociation from cohesin.
So a knob in this system can change motor state, not only add loss.

"Measured" means the design records the item directly. No such dataset has been obtained or checked.

| Contract item | What this assay records | Status | Risk / what is needed |
|---|---|---|---|
| Observed quantity | Presence of a labeled cohesin spot on each tether at one time | Measured | The "loss" event must mean cohesin leaving DNA. A knob that only stops extrusion (Wutz 2026, PDS5) produces bound, stalled complexes, not loss. |
| Independent unit | Chamber × protein prep × day; tethers cluster within chambers | Inferable with assumptions | The repo's binomial intervals assume independent tethers. Needs chamber-level or hierarchical intervals and ≥6 chambers per arm. |
| Initiated denominator | Tethers exposed to the loading pulse, counted from DNA imaging | Measured | Multiple loadings per tether: define presence as ≥1 complex, or keep occupancy low (Poisson). |
| Retention (absolute) | Present / initiated tethers per arm | Measured | Detection sensitivity s is absorbed into ℓ; needs identical imaging across arms and target. |
| Loading reporter (t=0) | Reference chambers imaged right after the wash, before any knob | Measured | Sequential design makes κ=1 between arms up to chamber variation. The reference then bounds chamber variation instead of a treatment effect. Requires randomized chamber allocation (exchangeability). |
| Sensitivity / specificity | Sensitivity = labeling efficiency × detection; false positives = stuck dye | Sensitivity: inferable with assumptions (bleaching-step counts, dual-label colocalization). False positives: measured (no-NIPBL / no-ATP controls). | q1/q0=1 by design, because the reference precedes the knob. A false-positive rate ≤1% has negligible impact here (§1). |
| Perturbation factor f | Dose ratio set by pipetting | Inferable with assumptions | Multiplier = dose ratio only if the hazard is linear in dose. The zero-dose arm measures δ alone; doses 0/1/2 test linearity. Needs ±1–2% (§2). |
| Target multiplier g | Dose 3 | **Missing** | The theory assumes g exact. In practice it has the same linearity assumption as f, extrapolated one step further, and the target cannot be used to test it. |
| Observation clock | Dark interval per chamber, Exp(γ) or scheduled | Measured | One frame per field keeps photobleaching before observation near zero, so δ is dose-independent dissociation, measured by the zero-dose arm. |
| Equal susceptibility (h=1) | Not visible in presence-only data | **Missing** | Extruding, stalled, CTCF-blocked and entrapped complexes may be released at different rates. Wutz 2026 shows a state-specific PDS5 effect. Needs state-resolved dwell times. h≠1 breaks the reciprocal law; the repo already demonstrates this. |
| Direct dwell-time readout | Interval imaging of the same tethers | Measured (natural in this assay) | The strong control. It distinguishes the 9/286 worlds (calibration-arm hazards 3/5 vs 7/10) and needs neither κ nor q. Two frame intervals separate bleaching (Reisser et al. 2020 precedent). |
| Held-out target cohort | Dose-3 chambers, labels blinded until intervals are locked | Measured | Allocate before any calibration analysis. |

**What the table shows.** This assay can meet the contract, but by design rather than through the
counted reference. Sequential loading removes κ, and a zero-dose arm removes the interval on f. The
target becomes a Stern–Volmer dose extrapolation, and dwell-time imaging dominates presence-only
snapshots. The assays that would *need* the counted reference are those where loading co-varies with
the perturbation and only one destructive presence readout per episode is possible. Live-cell assays
fit that description but lack initiated denominators, as the review packet already records. No
eligible assay with that combination has been identified.

## 4. Chosen next experiment: reconstitution feasibility pilot (no target data)

This is a design for a collaborator, frozen before any data. It touches only doses 0, 1 and 2. Dose 3
stays untouched.

Arms and allocation: doses 0/1/2, with ≥6 chambers per arm across ≥2 protein preps. In every chamber
set: t=0 reference chambers, single-frame snapshot chambers at randomized dark times, and
interval-imaging chambers at two frame intervals. Also run no-NIPBL and no-ATP specificity controls.
Allocate chambers at random, before imaging.

Go criteria (all required):
- **G1, the knob removes cohesin:** the 95% CI for the dose-2 minus dose-0 hazard excludes 0. The
  regulated share at dose 1, μ/(μ+δ), is ≥0.25; otherwise the target is insensitive to dose.
- **G2, linearity:** the 95% CI for the realized multiplier (a2−a0)/(a1−a0) lies within [1.96, 2.04]
  (±2%).
- **G3, equal susceptibility:** at every dose, the hazards of extruding and non-extruding complexes
  agree within ±10%, and the long-lived fraction does not change with dose. A two-component dwell
  model must not beat one component on held-out chambers by more than a margin set in advance.
- **G4, loading exchangeability:** the 95% CI for the between-arm loading ratio from the t=0 references
  lies within [0.98, 1.02].
- **G5, specificity:** false-positive spots on control tethers are ≤1% of tethers.

Stop rules: if G3 fails, the h=1 target does not apply to this knob; record it as a negative
applicability result and keep the paper theoretical. If G2 fails, a dose-defined target needs a
nonlinear model, which is outside the current theory. If everything passes, freeze the
dose-3 prediction from dose 0–2 data alone, using both designs: snapshots with the usefulness rule
from §2, and dwell times. Then unblind dose 3.

Why not a cost-matched simulation first: under this assay's contract, the exact-data analysis already
predicts the main outcome. With sequential loading and a zero-dose arm, structural width is zero
under linearity, so the reference adds nothing. The open questions (G2, G3, the definition of loss)
are biological. A simulation is worth running only for a declared assay where κ≠1 cannot be avoided,
and none is in hand.

## 5. Paper consequence

Keep the contribution narrow and theoretical: when a loading reference is *necessary* (loading
co-varies with the perturbation, readout is presence-only and destructive); the sharp drift
polytope and its corrected scope; and the calibration tolerances in §2. The reconstitution example
shows when the reference is *unnecessary*. Do not claim a protocol advantage over dwell-time imaging
or zero-dose Stern–Volmer designs. External proof and priority review is still required.

## Reproduce

```sh
python audit_counted_contract.py > results/counted_contract_audit.json
python -m unittest discover -s tests -p test_counted_contract_audit.py -v
```

[`counted_contract_audit.json`](../results/counted_contract_audit.json) holds every number quoted
above. `audit_counted_contract.py` imports nothing from the repository, uses no optimizer, and uses
exact `fractions` arithmetic for every reported check.
