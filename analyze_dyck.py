"""Where do models fail on Dyck? Accuracy by the depth of material to skip.

For a position whose innermost open bracket is at o, "skip depth" is the deepest
nesting among the completed pairs between o and the position (0 = only filler
in between). A model that collapses innermost loops one layer at a time should
degrade once skip depth exceeds what its layers can peel away.

    python train.py --task dyck --mixer extrusion --seed 0 --save results/checkpoints/dyck_extrusion_s0.pt
    python analyze_dyck.py extrusion one_arm bilstm attn_alibi
    python analyze_dyck.py extrusion_6L one_arm_6L   # checkpoints saved as dyck_<name>_s0.pt
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch

import task_dyck as T
from models import SequenceModel
from train import TEST_LENGTHS


def skip_depths(seq):
    """Per position: -1 if no bracket is open, else the skip depth."""
    out, stack = [], []  # stack entries: [type, deepest completed pair inside so far]
    for t in seq:
        if t in T.CLOSER:
            stack.append([t, 0])
        elif t != T.DOT:
            _, inner = stack.pop()
            if stack:
                stack[-1][1] = max(stack[-1][1], inner + 1)
        out.append(stack[-1][1] if stack else -1)
    return np.array(out)


@torch.no_grad()
def by_depth(model):
    hits, totals = defaultdict(int), defaultdict(int)
    for n, (seed, count) in TEST_LENGTHS.items():
        x, y = T.make_split(seed, count, [n])[n]
        pred = torch.cat([model(torch.from_numpy(x[s:s + 16])).argmax(-1) for s in range(0, len(x), 16)]).numpy()
        for seq, lab, p in zip(x, y, pred):
            for d, ok in zip(skip_depths(seq), p == lab):
                hits[int(d)] += int(ok)
                totals[int(d)] += 1
    return {d: (hits[d] / totals[d], totals[d]) for d in sorted(totals)}


def main(mixers, out_path="results/dyck_skip_depth.json"):
    torch.set_num_threads(4)
    table = {}
    for mixer in mixers:
        state = torch.load(f"results/checkpoints/dyck_{mixer}_s0.pt")
        layers = 1 + max(int(k.split(".")[1]) for k in state if k.startswith("blocks."))
        model = SequenceModel(mixer.removesuffix(f"_{layers}L"), vocab=len(T.TOKENS), n_classes=T.N_CLASSES, layers=layers)
        model.load_state_dict(state)
        model.eval()
        table[mixer] = by_depth(model)
    depths = sorted({d for rows in table.values() for d in rows})
    print("skip depth | positions | " + " | ".join(mixers))
    for d in depths:
        count = next(rows[d][1] for rows in table.values() if d in rows)
        label = "no open bracket" if d == -1 else str(d)
        print(f"{label:>10} | {count:>9} | " + " | ".join(f"{table[m][d][0]:.3f}" if d in table[m] else "-" for m in mixers))
    out = {m: {str(d): {"accuracy": round(a, 4), "positions": c} for d, (a, c) in rows.items()} for m, rows in table.items()}
    Path(out_path).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    names = sys.argv[1:] or ["extrusion", "one_arm", "bilstm", "attn_alibi"]
    deep = any(n.endswith("L") for n in names)
    main(names, "results/dyck_skip_depth_6layer.json" if deep else "results/dyck_skip_depth.json")
