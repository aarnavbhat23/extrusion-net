# Scan execution and physical parameter tying

Protocol written before the new training runs. Existing generality results are
already known and remain unchanged.

## Implementation question

Does the complete unsaturated mixer use linear-size computation in both forward
and backward, while preserving dense outputs and parameter/input gradients?

The earlier standalone recurrence was not wired into the mixer. Its repeated
slice writes also caused avoidable full-output copies in automatic differentiation.
Use stacked recurrence states and scan a constant channel to obtain contact mass.
Preserve the dense contact routine as a small reference and the saturated path.

Tests compare every parameter gradient, exact hard barriers, singleton sequences,
saved activation growth, and execution with dense contact construction forbidden.
The scaling script measures the complete mixer, not just the isolated recurrence.
Report Python-loop overhead and do not imply an optimized device kernel.

## Scientific control

The generic `cascaded_scan` control has independently learned load and pass gates
for each triangle. Tying them exactly reproduces linear extrusion. The extra
parameters are reported. This tests whether the physical constraint helps compared
with a more general member of the same structured-matrix family. It is not an
official implementation of Hydra or Mamba-2 and must not be described as one.

## Fixed pilot

- Task: existing synthetic genome task, 300 training sequences, one layer.
- Mixers: `extrusion_linear`, `cascaded_scan`, `attn_alibi`.
- Seeds: 0, 1, 2. Training: 500 updates, batch 32, learning rate 0.002;
  remaining optimizer/data/validation/test settings are those in `train.py`.
- Select checkpoints by short-length validation only. Save every run and checkpoint.
- Report accuracy at all four test lengths, parameter counts, and elapsed time.
- This is a short-update pilot, not a rerun of the published 3000-update grid.
- A preliminary advantage requires the physical model to exceed both controls in
  mean length-512 accuracy, with the per-seed comparisons visible. Three seeds
  cannot establish statistical significance or broad superiority.
- If the generic cascade matches or exceeds extrusion, do not claim a measured
  advantage for the physical constraint in this setting. Do not tune from the
  result or create a new mechanism to rescue it within this pilot.

The task's labels follow extrusion-like rules. A positive result is evidence of
an inductive bias on that task, not biological validation or a new learning rule.

## Novelty scope

The scan is a correct implementation of the derivation already in this branch;
it does not establish a new matrix class. Loop traversal needs a separate audit
against differentiable stacks, structured parsing, and learned graph rewrites
before any novelty or context-free capability claim. A simulation of colliding
motors likewise requires independent validation, rather than changing the label
generator to exactly match the proposed architecture.
