"""Train every mixer on short sequences, then test on longer ones.

Checkpoint selection uses an in-distribution validation set only; the longer
test lengths are never seen during training or selection.

    python train.py --mixer extrusion --seed 0
"""

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

import task
from models import SequenceModel, n_params

TRAIN_LENGTHS = [24, 32, 40, 48, 56, 64]
TEST_LENGTHS = {64: (3, 300), 128: (4, 300), 256: (5, 200), 512: (6, 100)}


def to_tensors(split):
    return {n: (torch.from_numpy(x), torch.from_numpy(y)) for n, (x, y) in split.items()}


@torch.no_grad()
def accuracy(model, split, chunk=16):
    model.eval()
    correct = total = 0
    for x, y in split.values():
        for s in range(0, len(x), chunk):
            xb, yb = x[s:s + chunk], y[s:s + chunk]
            pred = model(xb).argmax(-1)
            keep = yb != task.IGNORE
            correct += (pred[keep] == yb[keep]).sum().item()
            total += keep.sum().item()
    model.train()
    return correct / total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mixer", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-train", type=int, default=4000)
    ap.add_argument("--layers", type=int, default=3)
    ap.add_argument("--steps", type=int, default=3000)
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--lr", type=float, default=2e-3)
    ap.add_argument("--out", default="results/runs.jsonl")
    ap.add_argument("--save", help="optional path for the selected checkpoint")
    args = ap.parse_args()

    torch.set_num_threads(2)
    train = to_tensors(task.make_split(1, args.n_train, TRAIN_LENGTHS))
    val = to_tensors(task.make_split(2, 600, TRAIN_LENGTHS))
    tests = {n: to_tensors(task.make_split(s, k, [n])) for n, (s, k) in TEST_LENGTHS.items()}

    torch.manual_seed(args.seed)
    rng = np.random.default_rng(100 + args.seed)
    model = SequenceModel(args.mixer, layers=args.layers)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.01)
    warmup = 100
    sched = torch.optim.lr_scheduler.LambdaLR(
        opt,
        lambda t: min(1.0, (t + 1) / warmup)
        * 0.5 * (1 + math.cos(math.pi * min(t, args.steps) / args.steps)),
    )

    lengths = sorted(train)
    weights = np.array([len(train[n][0]) for n in lengths], dtype=float)
    weights /= weights.sum()
    best_val, best_state, curve = -1.0, None, []
    start = time.time()
    for step in range(1, args.steps + 1):
        n = lengths[rng.choice(len(lengths), p=weights)]
        x, y = train[n]
        idx = torch.from_numpy(rng.choice(len(x), size=min(args.batch, len(x)), replace=False))
        logits = model(x[idx])
        loss = F.cross_entropy(logits.reshape(-1, 4), y[idx].reshape(-1), ignore_index=task.IGNORE)
        opt.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        sched.step()
        if step % 250 == 0 or step == args.steps:
            v = accuracy(model, val)
            curve.append((step, round(loss.item(), 4), round(v, 4)))
            if v > best_val:
                best_val = v
                best_state = {k: t.clone() for k, t in model.state_dict().items()}

    model.load_state_dict(best_state)
    if args.save:
        torch.save(best_state, args.save)
    result = {
        "mixer": args.mixer,
        "seed": args.seed,
        "n_train": args.n_train,
        "layers": args.layers,
        "steps": args.steps,
        "params": n_params(model),
        "val_acc": round(best_val, 4),
        "test_acc": {str(n): round(accuracy(model, s), 4) for n, s in tests.items()},
        "train_seconds": round(time.time() - start, 1),
        "curve": curve,
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "a") as f:
        f.write(json.dumps(result) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "curve"}))


if __name__ == "__main__":
    main()
