"""Aggregate results/runs.jsonl into results/RESULTS.md (mean and range over seeds)."""

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ORDER = ["extrusion", "extrusion_linear", "one_arm", "attn_rope", "attn_alibi", "bilstm"]


def main(path="results/runs.jsonl", out="results/RESULTS.md"):
    runs = [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]
    groups = defaultdict(list)
    for r in runs:
        groups[(r["mixer"], r.get("task", "genome"), r["n_train"], r.get("layers", 3))].append(r)
    lengths = sorted({int(k) for r in runs for k in r["test_acc"]})

    lines = ["# Results", "", "Promoter-label accuracy (4 classes). Mean over seeds, with min-max in brackets.", ""]
    tasks = ["genome", "dyck", "recall"]
    settings = sorted({k[1:] for k in groups}, key=lambda s: (tasks.index(s[0]), -s[2], -s[1]))
    for task_name, n_train, layers in settings:
        lines += [f"## {task_name}: {n_train} training sequences (lengths 24-64), {layers} layer{'s' * (layers > 1)}", ""]
        header = "| Mixer | Params | Seeds | " + " | ".join(f"len {n}" for n in lengths) + " |"
        lines += [header, "|" + "---|" * (3 + len(lengths))]
        for mixer in ORDER:
            rs = groups.get((mixer, task_name, n_train, layers))
            if not rs:
                continue
            cells = []
            for n in lengths:
                acc = np.array([r["test_acc"][str(n)] for r in rs])
                cells.append(f"{acc.mean():.3f} [{acc.min():.3f}-{acc.max():.3f}]")
            lines.append(f"| `{mixer}` | {rs[0]['params']:,} | {len(rs)} | " + " | ".join(cells) + " |")
        lines.append("")
    Path(out).write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main(*sys.argv[1:])
