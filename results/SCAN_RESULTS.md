# Scan execution and physical tying: results

The study protocol is in `docs/SCAN_STUDY.md`. It was written locally before these runs; the code/protocol commit was made with the pilot already complete. This was not an externally preregistered trial. The polynomial-study protocol was committed before that numerical study. No test-selected hyperparameter changes were made.

## Fixed 500-update pilot

The existing genome task, 300 training sequences and one layer were used. All runs use the same data, optimizer, validation selection and update count. Means and ranges use three seeds; no significance claim is made.

| Mixer | Parameters | Length 64 | Length 128 | Length 256 | Length 512 | Individual seed accuracies at 512 |
|---|---:|---:|---:|---:|---:|---|
| extrusion_linear | 47,700 | 0.9497 | 0.9560 | 0.9515 | 0.9526 | 0.9050, 0.9767, 0.9762 |
| cascaded_scan | 48,480 | 0.9335 | 0.9336 | 0.9321 | 0.9284 | 0.9339, 0.9189, 0.9325 |
| attn_alibi | 54,980 | 0.7480 | 0.7201 | 0.6936 | 0.7017 | 0.6877, 0.6991, 0.7184 |

The physically tied model has a higher mean than both controls, satisfying the stated preliminary mean criterion. It loses to the untied cascade on seed 0. The cascade contains the tied model exactly and has 1.6% more parameters. The task follows an extrusion-like rule, so this supports a narrow inductive-bias claim. It does not establish superiority to official Hydra/Mamba implementations, biological performance, or a new learning paradigm.

All nine checkpoints reproduce their recorded accuracies at all four test lengths (`results/scan_pilot/replay_audit.json`). Token-level predictions and targets are stored with checksums. The original saturated checkpoint's replay also passes with unchanged accuracies and motor orientation.

## Complete-mixer execution

Seventeen tests pass. In particular, every input/parameter gradient matches a dense reference; exact hard barriers and singleton sequences are checked; the full linear model executes when dense contact construction is forbidden; and saved activation growth is linear over the checked sequence lengths.

The complete-mixer benchmark uses one CPU thread, width 32, two heads and batch one. Three timed repetitions follow a warm-up. Timing is machine-specific and includes Python-loop overhead. Logical saved-tensor bytes count views repeatedly and are not peak resident memory.

| Length | Path | Forward seconds | Forward + backward seconds | Logical saved bytes |
|---:|---|---:|---:|---:|
| 128 | scan | 0.00216 | 0.01286 | 213,280 |
| 128 | dense | 0.00122 | 0.00339 | 1,049,616 |
| 512 | scan | 0.00878 | 0.05637 | 827,680 |
| 512 | dense | 0.01838 | 0.04728 | 15,574,032 |
| 1024 | scan | 0.01676 | 0.11753 | 1,646,880 |
| 1024 | dense | 0.07851 | 0.19907 | 61,547,536 |
| 4096 | scan | 0.06790 | 0.64429 | 6,562,080 |

At short lengths the dense implementation is faster, particularly in backward. At 1024 the scan is faster in both measurements and stores far fewer logical activations. The benchmark supports the linear execution correction, not a claim about a fused GPU kernel or universal speedup. Dense runs above 1024 were deliberately omitted from this complete-gradient benchmark to limit memory use.

## Bounded saturation experiment

Exact entrywise contact-power scans and their gradients match dense references. In the seeded length-256 weak-contact example, order eight gives a maximum absolute message error of about 5.0e-9, below the analytic exact-arithmetic bound of 1.2e-5. Its Python implementation is slower than the dense reference at that length; the recorded one-shot polynomial timings are descriptive and include warm-up effects, not a speed benchmark.

The all-open length-32 example at beta=4 has maximum message error about 5.1e12 at order eight. A probe of the released checkpoint's first mixer has actual arguments up to 12.93 and much larger conservative bounds. The demonstrated weak-argument regime therefore does not justify enabling this approximation in that checkpoint. These failures are part of the result.

The specialized recurrence may be useful after a numerically stable approximation or a defensible bounded parameterization is developed. Presently it is a tested mathematical prototype with substantial prior overlap, not an improved trained architecture.

A separate inference-only check of the complete unsaturated mixer at 65,536 positions, width 32, two heads and batch one took 0.853 seconds on one CPU thread with finite outputs. This is one timing, not full-network training (`results/long_mixer_execution.json`).
