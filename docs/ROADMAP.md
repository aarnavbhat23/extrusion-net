# Roadmap: from one nature-inspired operator to a model worth scaling

## The lens: does it collapse?

Many nature-inspired ideas look new but turn out to be a known algorithm in disguise once written as math. Each candidate here was screened by writing its computation down and asking what it reduces to:

| Nature system | What it computes | Reduces to | Verdict |
|---|---|---|---|
| Bacterial bet-hedging (Kussell and Leibler, Science 2005) | switch phenotypes at the rate the environment changes | fixed-share expert tracking (Herbster and Warmuth 1998) | collapses |
| Germinal centers, antibody epitope masking | solved epitopes are hidden, forcing new coverage | boosting's reweighting | collapses |
| Consumer-resource ecology (MacArthur 1970) | species abundances at equilibrium | a constrained quadratic minimization | collapses |
| Kinetic proofreading (Hopfield, PNAS 1974) | accuracy from several reset-prone steps | a sharpened nonlinearity / sequential test | collapses |
| Hsp90 cryptic variation (Rutherford and Lindquist, Nature 1998) | hidden variants released under stress | ensembles with selection after a shift | collapses |
| Whi3 assemblies (Caudron and Barral, Cell 2013) | persistent memory from assemblies | a memory gate (found earlier in this project) | collapses in its simple form |
| Transcription supercoiling (Liu and Wang, PNAS 1987) | signed, directional twist spreading to barriers | signed bidirectional state-space scan | collapses |
| Morphogen gradients | source, diffusion, decay | convolution / state-space scan | collapses |
| **Loop extrusion, one motor** (this repo) | contact via a load site and two oriented arms | linear part is 2-quasiseparable; saturation and parameterization are not standard | **partly**: measurable advantages survived |
| **Loop traversal** (motors passing each other's loops; Kim et al., Nature 2020) | arms skip over completed loops | stack-like matching; no known operator family | **open** |
| **Condensates with a nucleation threshold** (phase separation; Shin and Brangwynne, Science 2017) | content-based grouping that forms only above a quorum, with competition for a limited scaffold | slot attention, plus a threshold that may be new | **open, high collapse risk** |

The candidates still standing come from the same biological system: how the nucleus decides which parts of the genome touch. Cells use two independent mechanisms for that: loop extrusion and compartmental segregation (Schwarzer et al., Nature 2017; Nuebler et al., PNAS 2018). That points to a family of operators rather than a single one.

## Candidate directions, ranked

### A1. Extrusion at genome scale on real enhancer-gene data (highest priority)

The linear-time scan (`extrusion_linear_scan`) handles 65,536 positions on one CPU thread. That makes the real question testable: does extrusion predict which enhancers regulate which genes better than distance and activity-by-contact (ABC) baselines? Data: CRISPRi enhancer-gene screens collected in the ENCODE-rE2G benchmark (Fulco et al., Nature Genetics 2019; Gschwind et al., bioRxiv 2023). Inputs: binned 1-2 Mb windows with CTCF motif orientation and strength, enhancer activity (DNase, H3K27ac) and gene TSS positions. Split by chromosome.

A win here is a biology result, not a toy result, and it is the step most likely to make the project matter.

### A2. Loop traversal for nested structure

The current model has one motor at a time. In cells, motors meet: they block each other, or one passes another and forms a Z-loop. If an arm can step over a loop that is already complete, it skips matched pairs, which is exactly what a stack does. That could give context-free capability (brackets, syntax, RNA secondary structure) from a physical rule rather than a learned stack. A differentiable version would be iterative: each layer forms the innermost loops, then marks their anchors as passable for the next layer, so nesting depth costs layers rather than memory.

Motivated by the Dyck result: extrusion reached 0.987 at length 512, above attention but 0.4 points below BiLSTM, and every model degraded on deeply nested material. Doubling layers did not help, so stacking is not a substitute for a skip rule. The decisive test is depth generalization: train with nesting depth at most 3, test at depths 4-8. A physical rule that steps over completed loops should not care about depth; learned stacks usually do. RNA secondary structure, which is also nested, is a natural real-data follow-up and connects to `rna-shift-bench`.

### A3. A "nuclear computer": extrusion plus condensates

Extrusion has no lookup by content; it cannot match a query to an arbitrary key. The genome's other contact mechanism can: phase-separated condensates gather elements by shared binding affinity. A model with both mixers, extrusion for topology and condensates for content, would mirror the nucleus. The risk is that the condensate half collapses to attention. The testable claim is the combination: tasks that need both topology (insulation) and content (specificity), such as enhancer-promoter compatibility classes inside insulated neighborhoods (Zabidi et al., Nature 2015).

Motivated by the recall result: extrusion falls toward chance at length 512 while ALiBi attention stays at 0.99.

## Validation ladder (applies to every candidate)

Each stage has a pass condition decided before running it. Failing a stage means reporting it and stopping, not tuning until it passes.

| Stage | Question | Pass condition | Extrusion status |
|---|---|---|---|
| 0. Mechanism | Is the math correct, and is the bio rule representable? | exactness tests pass | passed |
| 1. Matched synthetic | Does the prior help where it should? | beats generic mixers at matched size | passed: length, data and step efficiency |
| 2. Collapse | Does it reduce to a known operator? | an advantage over the nearest known operator survives | partly: linear part is quasiseparable; beats one-arm (SSM-like) at 1 layer and at low data |
| 3. Generality | Where does it fail? | failure modes documented and explained | done: strong on nested brackets (0.4 points below BiLSTM), fails content lookup |
| 4. Mechanism mismatch | Does it survive noisy, stochastic data? | advantage holds on simulated multi-motor data | not started |
| 5. Real data | Does it beat domain baselines? | beats distance and ABC on held-out chromosomes | not started (A1) |
| 6. Scale | Does it run at real sizes? | linear-time version keeps the advantage | passed on toy data: within 0.9 points, 65,536 positions in 2 s |
| 7. Independent check | Does someone else agree? | external critique and replication | not started |

## Where outside help fits

- **Independent technical review:** useful at stage 2 and 7. Give it the derivation and ask it to find prior work on chained bidirectional scans and traversal-based parsing, and to attack the claims.
- **A mentor or lab:** needed at stage 5, to check the enhancer-gene data processing, and later for any wet-lab validation.

## Not claimed

None of this is a route to general intelligence, and "the next neural network" is a goal, not a result. The defensible claim today: a chromatin-inspired routing operator with measured advantages on tasks with extrusion-like structure, a linear-time form, and a documented relationship to known operators.
