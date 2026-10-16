# Extrusion Net

[![Reproduction checks](https://github.com/aarnavbhat23/extrusion-net/actions/workflows/reproduce.yml/badge.svg)](https://github.com/aarnavbhat23/extrusion-net/actions/workflows/reproduce.yml)

Can a neural network route information the way chromatin does?

In the nucleus, cohesin motors load onto DNA and reel it into loops from both sides. Each side stalls at a CTCF protein facing the right way, so a forward site (`>`) followed by a reverse site (`<`) anchors a loop that insulates genes from enhancers outside it. This repository turns that process into a sequence model's only long-range operator. Every position predicts, from its local context alone, how likely a motor is to load there and how likely each arm is to stall there. Two positions exchange information only if a simulated motor brings them into the same loop. There is no positional encoding and no query-key product.

On a toy genome task that follows the convergent-CTCF rule, the extrusion model trained on 4,000 sequences **keeps 99.8-100% accuracy at 8x its longest training sequence**. RoPE attention falls to 58%. With **300 training sequences** it reaches **98.4-99.0%**, while BiLSTM and the one-arm ablation reach 93.7-95.3% and attention 53-79%. Without being told which tokens are barriers, it **learned the CTCF orientation rule**.

![Accuracy by test length](results/generalization.png)

## What happened

36 runs: five long-range mixers on one shared backbone, 3 seeds each. Promoter accuracy, mean of 3 seeds; all models trained on lengths 24-64 only.

| Mixer | Params | 4,000 seqs, len 64 | 4,000 seqs, len 512 | 300 seqs, len 64 | 300 seqs, len 512 |
|---|---:|---:|---:|---:|---:|
| Extrusion (two-arm) | 141,428 | **1.000** | **1.000** | **0.984** | **0.986** |
| One-arm ablation | 141,428 | 0.998 | 0.997 | 0.953 | 0.937 |
| Attention + RoPE | 163,268 | 0.981 | 0.582 | 0.787 | 0.533 |
| Attention + ALiBi | 163,268 | 0.952 | 0.928 | 0.789 | 0.732 |
| BiLSTM | 188,612 | 0.997 | **1.000** | 0.950 | 0.940 |

With **one layer** (4,000 sequences), two-arm extrusion scores 0.998 at length 512 and the one-arm ablation 0.968. It wins on all 3 seeds at every length. All lengths, seed ranges and the 1-layer table are in [`results/RESULTS.md`](results/RESULTS.md).

1. **Length generalization.** Extrusion and BiLSTM stay near-perfect at 512. RoPE attention collapses; ALiBi plateaus around 93%.
2. **Sample efficiency.** With 300 sequences, at every test length, extrusion's worst seed beats every other model's best seed.
3. **The two-arm mechanism matters.** One extrusion layer can represent the convergent rule exactly ([proof](docs/DERIVATION.md#exactness-for-the-convergent-ctcf-rule)); one one-arm layer cannot. With 3 layers the ablation catches up through depth.
4. **With enough data, it ties BiLSTM.** The rule can be checked by a left-to-right scan, which suits recurrent models. Extrusion's advantage over BiLSTM here is data efficiency and parallel evaluation, not peak accuracy.

## What the model learned

![Learned motor probabilities](results/motors.png)

The model only ever saw promoter labels. In the released checkpoint, three of four first-layer motors learned that `>` stalls left-moving arms (0.995-0.999) but not right-moving ones (0.019-0.068). They also learned that `<` stalls right-moving arms (0.985-0.999) but not left-moving ones (0.005-0.023). That is the orientation rule found experimentally in real genomes. Loading is high at enhancers in all four motors (0.58-0.96), and in motors 2 and 4 also at forward sites (0.96-0.99). Motor 3 learned something else. Only this one checkpoint was inspected.

## The mechanism

For positions `i < j`, with load propensity `l`, left-arm stall `sL` and right-arm stall `sR`:

```
D(i, j) = sum_{k=i+1..j} l_k * prod_{m=i+1..k-1} (1 - sL_m) * prod_{m=k..j-1} (1 - sR_m)
W(i, j) = 1 - exp(-beta * D(i, j))
```

Position `i` receives `sum_j W(i, j) v_j`. Prefix sums make this exact in O(n^2), the same order as attention. The derivation, numerics and exactness proof are in [`docs/DERIVATION.md`](docs/DERIVATION.md).

## Reproduce

Python 3.13 or 3.14 on a CPU. The recorded runs used an Apple M3, 2 threads per run.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python replay.py
python predict.py "3..P..>..<..2..P"
```

`replay.py` reloads the released checkpoint and checks its recorded test accuracy and learned orientation motors. `predict.py` takes any string of `. 1 2 3 P > <` and reports each promoter's prediction next to the rule. To rerun experiments:

```bash
python train.py --mixer extrusion --seed 0                         # one row of the grid, about 10 CPU minutes
python train.py --mixer one_arm --seed 0 --layers 1                # 1-layer setting
python train.py --mixer bilstm --seed 0 --n-train 300              # low-data setting
python summarize.py && python inspect_learned.py && python plot.py
```

Mixers: `extrusion`, `one_arm`, `attn_rope`, `attn_alibi`, `bilstm`. Every run appends one line to `results/runs.jsonl`.

## Protocol

Task, splits, seeds, models, training and selection are specified in [`docs/PROTOCOL.md`](docs/PROTOCOL.md). In short: tokens `. E1 E2 E3 P > <` are sampled independently. Each promoter's label is the strongest enhancer (0-3) not insulated from it by a `>` preceding a `<`. Train lengths are 24-64; test lengths are 64, 128, 256 and 512. Checkpoints are selected on in-distribution validation only. All five mixers share the embedding, convolution, MLP and head.

That document also says which choices were made after seeing results. The 1-layer and 300-sequence settings were added because the first grid saturated. The released checkpoint was retrained for inspection, so its numbers differ slightly from the grid's seed-0 row.

## Artifacts

- `results/runs.jsonl`: all 36 runs, with validation curves.
- `results/RESULTS.md`: per-setting tables with seed ranges.
- `results/extrusion_seed0.pt`: released extrusion checkpoint.
- `results/replay.json`: its recorded test accuracy and layer-1 orientation motors.
- `results/motors.json`, `results/inspect.txt`: learned stall and load probabilities by token, all layers.
- `results/*.png`, `results/*.svg`: figures from `plot.py`.

## Prior work and credit

Loop extrusion and the convergent-CTCF rule come from chromatin biology (Sanborn et al. 2015; Fudenberg et al. 2016; Rao et al. 2014; de Wit et al. 2015; Davidson et al. 2019; Kim et al. 2019). The closest machine-learning work is dLEM (2025), which makes extrusion differentiable to predict Hi-C contact maps, and stick-breaking attention (Tan et al., ICLR 2025), whose cumulative stopping resembles the one-arm ablation. Full comparison, and the searches not yet done, are in [`docs/PRIOR_WORK.md`](docs/PRIOR_WORK.md). No first-of-its-kind claim is made.

## Limitations

- **This is a toy task built on the mechanism's own rule**, like testing a convolution on translated images. It shows that the prior helps when data follow it, not that the model is better in general.
- Labels are deterministic and symbolic. Real chromatin has leaky barriers, finite processivity, motor collisions and noisy measurements.
- 3 seeds per setting; no significance tests. One checkpoint was inspected.
- O(n^2) memory, like attention. Genome-scale inputs would need a banded version based on limited motor travel.

## Next experiments

1. Labels generated by a stochastic multi-motor extrusion simulation, so the model is an approximation rather than an exact match.
2. Tasks outside genomics with nested structure (Dyck languages, ListOps) under the same length protocol.
3. Real data: enhancer-gene links from CRISPRi screens (the ENCODE-rE2G benchmark), using CTCF motif orientation, against distance and activity-by-contact baselines.
4. Ablations: shared stall probability for both arms (no orientation), no Poisson saturation, no contact-mass channel.

Stop rules for these are fixed in [`docs/PROTOCOL.md`](docs/PROTOCOL.md#7-fixed-before-the-next-experiments).

Built with AI assistance. All numbers above come from the recorded runs, not generated example outputs. Contributions that reproduce results or identify a failure are welcome; there is no independent review or adoption claim.
