"""Did the model discover the CTCF orientation rule from labels alone?

Averages each motor type's learned stall and load probabilities by token for
the released extrusion checkpoint. The model only ever saw promoter labels; it
was never told which token is a barrier or which direction a barrier faces.

    python inspect_learned.py           # uses results/extrusion_seed0.pt
    python inspect_learned.py --retrain # retrains it first (about 10 CPU minutes)
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

import task
from models import SequenceModel
from train import TRAIN_LENGTHS, accuracy, to_tensors

CHECKPOINT = Path("results/extrusion_seed0.pt")


def train_model(steps=3000, seed=0):
    """Same data, schedule and selection rule as train.py --mixer extrusion --seed 0."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(100 + seed)
    train = to_tensors(task.make_split(1, 4000, TRAIN_LENGTHS))
    val = to_tensors(task.make_split(2, 600, TRAIN_LENGTHS))
    model = SequenceModel("extrusion")
    opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=0.01)
    sched = torch.optim.lr_scheduler.LambdaLR(
        opt, lambda t: min(1.0, (t + 1) / 100) * 0.5 * (1 + math.cos(math.pi * min(t, steps) / steps))
    )
    lengths = sorted(train)
    weights = np.array([len(train[n][0]) for n in lengths], dtype=float)
    weights /= weights.sum()
    best, best_state = -1.0, None
    for step in range(1, steps + 1):
        x, y = train[lengths[rng.choice(len(lengths), p=weights)]]
        idx = torch.from_numpy(rng.choice(len(x), size=32, replace=False))
        loss = F.cross_entropy(model(x[idx]).reshape(-1, 4), y[idx].reshape(-1), ignore_index=task.IGNORE)
        opt.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        sched.step()
        if step % 250 == 0:
            v = accuracy(model, val)
            if v > best:
                best, best_state = v, {k: t.clone() for k, t in model.state_dict().items()}
    model.load_state_dict(best_state)
    return model


@torch.no_grad()
def motor_tables(model, n_seqs=300, n=64):
    """Per layer: {token: [[stall_left, stall_right, load] per motor]}."""
    rng = np.random.default_rng(7)
    x = torch.from_numpy(np.stack([task.sample_sequence(rng, n) for _ in range(n_seqs)]))
    model.eval()
    h = model.embed(x)
    tables = []
    for block in model.blocks:
        h = h + F.gelu(block.conv(block.ln_conv(h).transpose(1, 2)).transpose(1, 2))
        mixer = block.mixer
        z = mixer.control(block.ln_mix(h)).view(n_seqs, n, mixer.heads, 3)
        probs = torch.sigmoid(z.clamp(max=12))
        tables.append({name: probs[x == t].mean(0).tolist() for t, name in enumerate(task.TOKENS)})
        h = h + mixer(block.ln_mix(h))
        h = h + block.mlp(block.ln_mlp(h))
    return tables


def format_tables(tables):
    lines = []
    for layer, table in enumerate(tables, 1):
        lines.append(f"Layer {layer}: mean learned probability by token (stall-left / stall-right / load)")
        heads = len(next(iter(table.values())))
        lines.append("token " + "".join(f"| motor {h + 1}            " for h in range(heads)))
        for name, p in table.items():
            lines.append(f"{name:5} " + "".join(f"| {a:.2f} / {b:.2f} / {c:.2f} " for a, b, c in p))
        lines.append("")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--retrain", action="store_true")
    args = ap.parse_args()
    torch.set_num_threads(4)
    if args.retrain or not CHECKPOINT.exists():
        model = train_model()
        torch.save(model.state_dict(), CHECKPOINT)
    else:
        model = SequenceModel("extrusion")
        model.load_state_dict(torch.load(CHECKPOINT))
    tables = motor_tables(model)
    Path("results/motors.json").write_text(json.dumps(tables, indent=1))
    text = format_tables(tables)
    Path("results/inspect.txt").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
