# Roadmap: from one nature-inspired operator to a model worth scaling

## Current selection: learn transport through observable probes

The [loading-probe screen](LOAD_PROBES.md) is now the primary learning-research
direction. The follow-up [global recovery theorem](COLORED_RECOVERY.md) now gives
a length-independent sufficient response budget for every fixed gate ceiling
below one, with a projected recovery algorithm and a posteriori noise bound.
It retains an explicit source-invariance counterexample and held-out validation
test. Next reduce the transparency-dependent response cost using local barrier
conditions and test mechanism transfer beyond the matched simulator.

The earlier screen established that two controlled source patterns remove the free-gate ambiguity;
four vector responses offer a promising O(n)-observation route to stable local
recovery. Shared-controller fits transfer to new feature contexts, with random
input design a strong control. The new results settle a restricted operator
recovery question; a useful general learning advantage remains an open target.

The biological prediction track below remains useful for external validation.
The older architecture candidates are retained as hypotheses, not selected
implementations. Their analogies alone do not justify development.

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
| **Loop traversal** (motors passing each other's loops; Kim et al., Nature 2020) | arms skip over completed loops | stack-like matching and structured parsing; existing differentiable stacks are close | **open, high overlap risk** |
| **Condensates with a nucleation threshold** (phase separation; Shin and Brangwynne, Science 2017) | content-based grouping that forms only above a quorum, with competition for a limited scaffold | slot attention, plus a threshold that may be new | **open, high collapse risk** |

The candidates still standing come from the same biological system: how the nucleus decides which parts of the genome touch. Cells use two independent mechanisms for that: loop extrusion and compartmental segregation (Schwarzer et al., Nature 2017; Nuebler et al., PNAS 2018). That points to a family of operators rather than a single one.

## Candidate directions, ranked

### A1. Biological validation: real enhancer-gene data

The linear-time scan (`extrusion_linear_scan`) handles 65,536 positions on one CPU thread. That makes the real question testable: does extrusion predict which enhancers regulate which genes better than distance and activity-by-contact (ABC) baselines? Data: CRISPRi enhancer-gene screens collected in the ENCODE-rE2G benchmark (Fulco et al., Nature Genetics 2019; Gschwind et al., bioRxiv 2023). Inputs: binned 1-2 Mb windows with CTCF motif orientation and strength, enhancer activity (DNase, H3K27ac) and gene TSS positions. Split by chromosome.

A win here is a biology result, not a toy result, and it is the step most likely to make the project matter.

### A2. Loop traversal for nested structure

The current model has one motor at a time. In cells, motors meet: they block each other, or one passes another and forms a Z-loop. If an arm can step over a loop that is already complete, it skips matched pairs, which is exactly what a stack does. That could give context-free capability (brackets, syntax, RNA secondary structure) from a physical rule rather than a learned stack. A differentiable version would be iterative: each layer forms the innermost loops, then marks their anchors as passable for the next layer, so nesting depth costs layers rather than memory.

Motivated by the Dyck result: extrusion reached 0.987 at length 512, above attention but 0.4 points below BiLSTM, and every model degraded on deeply nested material. Doubling layers did not help, so stacking is not a substitute for a skip rule. The decisive test is depth generalization: train with nesting depth at most 3, test at depths 4-8. Depth extrapolation must be tested against established differentiable stacks and nondeterministic pushdown models; the physical analogy supplies no guarantee. RNA secondary structure, which is also nested, is a natural real-data follow-up and connects to `rna-shift-bench`.

### A3. A "nuclear computer": extrusion plus condensates

Extrusion has no lookup by content; it cannot match a query to an arbitrary key. The genome's other contact mechanism can: phase-separated condensates gather elements by shared binding affinity. A model with both mixers, extrusion for topology and condensates for content, would mirror the nucleus. The risk is that the condensate half collapses to attention. The testable claim is the combination: tasks that need both topology (insulation) and content (specificity), such as enhancer-promoter compatibility classes inside insulated neighborhoods (Zabidi et al., Nature 2015).

Motivated by the recall result: extrusion falls toward chance at length 512 while ALiBi attention stays at 0.99.

## Validation ladder (applies to every candidate)

Each stage has a pass condition decided before running it. Failing a stage means reporting it and stopping, not tuning until it passes.

| Stage | Question | Pass condition | Extrusion status |
|---|---|---|---|
| 0. Mechanism | Is the math correct, and is the bio rule representable? | exactness tests pass | passed |
| 1. Matched synthetic | Does the prior help where it should? | beats generic mixers at matched size | passed: length, data and step efficiency |
| 2. Collapse | Does it reduce to a known operator? | an advantage over the nearest known operator survives | partly: quasiseparable; a short untied-cascade control is now tested, with one seed losing; full structured-model comparisons remain |
| 3. Generality | Where does it fail? | failure modes documented and explained | done: strong on nested brackets (0.4 points below BiLSTM), fails content lookup |
| 4. Mechanism mismatch | Does it survive noisy, stochastic data? | advantage holds on simulated multi-motor data | not started |
| 5. Real data | Does it beat domain baselines? | beats distance and ABC on held-out chromosomes | not started (A1) |
| 6. Scale | Does it run at real sizes? | linear-time version keeps the advantage | complete unsaturated forward/backward now verified; 65,536-position number was the standalone recurrence |
| 7. Independent check | Does someone else agree? | external critique and replication | not started |

## Where outside help fits

- **Independent technical review:** useful at stage 2 and 7. Give it the derivation and ask it to find prior work on chained bidirectional scans and traversal-based parsing, and to attack the claims.
- **A mentor or lab:** needed at stage 5, to check the enhancer-gene data processing, and later for any wet-lab validation.

## Not claimed

None of this is a route to general intelligence, and "the next neural network" is a goal, not a result. The defensible claim today: a chromatin-inspired routing operator with measured advantages on tasks with extrusion-like structure, a linear-time form, and a documented relationship to known operators.

## Current gate after the execution/control study

Do not treat a correct scan as architectural novelty. The short tied/untied pilot supports a narrow inductive-bias hypothesis, with seed variability. The saturation approximation has a severe large-argument limitation. Before a traversal mechanism is built, compare its formal rule against [differentiable stacks](https://arxiv.org/abs/1506.02516) and [nondeterministic stack models](https://arxiv.org/abs/2010.04674). A stack under a biological name does not establish a new learning rule. Next evidence should come from a matched structured-model comparison or an appropriate biological task with controls, not additional analogies.
