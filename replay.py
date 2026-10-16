"""Recompute the released checkpoint's results and compare with the recorded ones.

Checks that results/extrusion_seed0.pt still gives the recorded test accuracy
at every length, and that its learned motor table still shows the CTCF
orientation rule. Fails loudly if anything drifts.

    python replay.py          # check
    python replay.py --write  # record (only after an intentional rerun)
"""

import argparse
import json
from pathlib import Path

import torch

import task
from inspect_learned import CHECKPOINT, motor_tables
from models import SequenceModel
from train import TEST_LENGTHS, accuracy, to_tensors

RECORD = Path("results/replay.json")
TOLERANCE = 0.002


def compute():
    model = SequenceModel("extrusion")
    model.load_state_dict(torch.load(CHECKPOINT))
    tests = {n: to_tensors(task.make_split(s, k, [n])) for n, (s, k) in TEST_LENGTHS.items()}
    acc = {str(n): round(accuracy(model, split), 4) for n, split in tests.items()}
    layer1 = motor_tables(model)[0]
    return {"test_acc": acc, "layer1_motors": {t: [[round(v, 3) for v in m] for m in layer1[t]] for t in (">", "<")}}


def orientation_motors(motors):
    """Motors where '>' stalls left arms and '<' stalls right arms, but not the reverse."""
    fwd, rev = motors[">"], motors["<"]
    return [h for h in range(len(fwd)) if fwd[h][0] > 0.95 and fwd[h][1] < 0.1 and rev[h][1] > 0.95 and rev[h][0] < 0.1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()
    torch.set_num_threads(4)
    got = compute()
    if args.write:
        RECORD.write_text(json.dumps(got, indent=1))
        print(f"recorded {RECORD}")
        return
    want = json.loads(RECORD.read_text())
    for n, value in want["test_acc"].items():
        diff = abs(got["test_acc"][n] - value)
        print(f"length {n}: recorded {value:.4f}, recomputed {got['test_acc'][n]:.4f}")
        assert diff <= TOLERANCE, f"length {n} drifted by {diff:.4f}"
    found = orientation_motors(got["layer1_motors"])
    print(f"layer-1 motors following the CTCF orientation rule: {[h + 1 for h in found]}")
    assert found == orientation_motors(want["layer1_motors"]), "learned orientation motors changed"
    assert len(found) >= 1, "no motor follows the orientation rule"
    print("replay OK")


if __name__ == "__main__":
    main()
