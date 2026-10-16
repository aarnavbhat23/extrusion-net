"""Ask the released model which enhancer can switch on each promoter.

One character per position: '.' filler, '1'/'2'/'3' enhancers, 'P' promoter,
'>' forward CTCF site, '<' reverse CTCF site. Any length works.

    python predict.py "3..P..>..<..2..P"
"""

import sys

import numpy as np
import torch

import task
from inspect_learned import CHECKPOINT
from models import SequenceModel

CHARS = {".": task.DOT, "1": task.E1, "2": task.E2, "3": task.E3, "P": task.PROM, ">": task.FWD, "<": task.REV}


def main(genome):
    unknown = sorted(set(genome) - set(CHARS))
    if unknown:
        sys.exit(f"unknown characters {unknown}; use only {''.join(CHARS)}")
    seq = np.array([CHARS[c] for c in genome])
    if not (seq == task.PROM).any():
        sys.exit("no promoter 'P' in the sequence")
    model = SequenceModel("extrusion")
    model.load_state_dict(torch.load(CHECKPOINT))
    model.eval()
    with torch.no_grad():
        pred = model(torch.from_numpy(seq)[None]).argmax(-1)[0].numpy()
    truth = task.label_sequence(seq)
    print(genome)
    for p in np.flatnonzero(seq == task.PROM):
        verdict = "ok" if pred[p] == truth[p] else "WRONG"
        print(f"promoter at {p:4d}: predicted strongest enhancer {pred[p]}, rule says {truth[p]}  [{verdict}]")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "3..P..>..<..2..P")
