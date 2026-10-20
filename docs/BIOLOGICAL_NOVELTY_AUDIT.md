# Biological-modeling novelty: scope and the published comparator

**Learnable loop extrusion is already an established research direction.** It
cannot itself be the novelty claim for this project. The strongest present
candidate contribution is controlled-probe identification, explicit ambiguity
and falsifiable intervention prediction under a specified observation model.
Priority for those narrower results remains unestablished.

A targeted primary-source audit on 2026-10-20 checked:

- [Chikina lab's dLEM description](https://chikinalab.org/tools.html): fitting
  loop-extrusion dynamics to Hi-C, integration with one-dimensional genomic
  features and perturbation prediction.
- [Authors' PyTorch implementation](https://github.com/chikinalab/dLEMpytorch):
  a deep model predicting contact maps from sequence and epigenomic features,
  with supplied checkpoints and inference instructions.
- [Authors' JAX forward operator](https://github.com/chikinalab/dLEM/blob/3063ee7cc1a39f2114ba0ffabaafb342793a801e/dlem/model.py),
  pinned at `3063ee7cc1a39f2114ba0ffabaafb342793a801e`: a recursive contact-band
  operator whose numerator transports the preceding band in two directions,
  divided by exit and detachment rates; self-normalized rollout preserves
  gradients during deeper parameter fitting.
- The [paper](https://doi.org/10.1101/2025.10.17.682904) was located, but its
  full-text endpoint returned HTTP 403. This audit relies on authors' code and
  descriptions; it does not establish that a capability is absent from the full
  paper. At the initial audit no dLEM performance comparison had been run.

The two codebases have different readouts and gate semantics. Extrusion's
current controlled experiments use expected cross-arm reach, while dLEM models
contact-band steady-state dynamics and fits chromatin contact experiments.
Raw prediction errors cannot be compared across these observables. Neither a
new label for the recurrence nor a different implementation proves superiority.
Before a competitive biological-model claim, define a common observation model,
training information, held-out perturbations and resource accounting, then run
the published implementation on that matched task.

The external [event-engine study](EXTERNAL_ENGINE_RESULTS.md) is a restricted
validation of a reach-model bridge. It is not evidence of beating dLEM, predicting
experimental Hi-C, or replacing neural networks or transformers. The general
content-retrieval candidate [failed its gates](MOMENT_MEMORY_RESULTS.md).

## Subsequent matched-observable test

The [random-delay snapshot study](OCCUPANCY_RESULTS.md) now defines a controlled
observable for which the unchanged authors' netdlem2 forward has an exact
resolvent mapping. Our memoryless recurrence agrees numerically with it. The
hidden-capture fit's 1.93× primary intervention-error improvement fails its 2×
gate. This compares one pinned forward module with all effective step rates
fitted, not the complete deep prediction pipeline or all dLEM variants. The
repository also contains a stalling module, not benchmarked here. A competitive
claim on measured contact maps remains unestablished. The comparison strengthens
the prior-art boundary rather than removing it.

## Research decision

Prioritize a specific question: **which affordable measurements make an
intervention prediction identifiable, and when should the model reject its own
mechanism?** Develop joint parameter/family uncertainty and compare it against
point fitting and existing mechanistic models at identical observation cost.
Require an independently generated or measured perturbation outcome. If a
matched comparison does not reveal a new capability or useful tradeoff, narrow
or stop the claim rather than add architectural components.

A future general biological architecture needs evidence across independently
specified biological systems, not only more variants of this transport simulator.
None is established here. The current theorems are conditional mathematical
results; broad foundational novelty should not be inferred from them.
