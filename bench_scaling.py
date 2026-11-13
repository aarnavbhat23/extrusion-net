"""Time unsaturated extrusion mixing: full n x n matrix vs the O(n) scan.

The scan here is a plain Python loop (not a parallel scan kernel), so absolute
times are pessimistic; the point is how each method grows with length.

    python bench_scaling.py
"""

import json
import time
from pathlib import Path

import torch

from models import extrusion_linear_scan, extrusion_log_contacts

DENSE_MAX = 4096  # beyond this the n x n intermediates need several GB


def inputs(n, d=16, seed=0):
    g = torch.Generator().manual_seed(seed)
    load = torch.rand(n, generator=g) * 0.5 + 0.05
    stall_l = torch.rand(n, generator=g) * 0.3
    stall_r = torch.rand(n, generator=g) * 0.3
    return load, 1 - stall_l, 1 - stall_r, torch.randn(n, d, generator=g)


def dense(load, pl, pr, v):
    w = extrusion_log_contacts(load.log(), pl.log(), pr.log()).exp()
    return w @ v


def timed(fn, *args):
    start = time.perf_counter()
    out = fn(*args)
    return out, time.perf_counter() - start


def main():
    torch.set_num_threads(1)
    rows = []
    for n in [256, 1024, 4096, 16384, 65536]:
        args = inputs(n)
        scan_out, scan_s = timed(extrusion_linear_scan, *args)
        row = {"length": n, "scan_seconds": round(scan_s, 3)}
        if n <= DENSE_MAX:
            dense_out, dense_s = timed(dense, *args)
            row["dense_seconds"] = round(dense_s, 3)
            row["max_rel_diff"] = float((dense_out - scan_out).abs().max() / dense_out.abs().max())
        row["dense_matrix_mb"] = round(n * n * 4 / 2**20, 1)
        rows.append(row)
        print(row)
    Path("results/scaling.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
