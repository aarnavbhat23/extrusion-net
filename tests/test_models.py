import unittest

import numpy as np
import torch

import task
from models import MIXERS, SequenceModel, extrusion_log_contacts, insulation_log_contacts


def brute_extrusion(load, pass_l, pass_r):
    n = len(load)
    d = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(i + 1, j + 1):
                d[i, j] += load[k] * np.prod(pass_l[i + 1:k]) * np.prod(pass_r[k:j])
            d[j, i] = d[i, j]
    return d


def brute_label(seq):
    """Direct statement of the rule, O(n^3), for checking task.label_sequence."""
    y = np.full(len(seq), task.IGNORE)
    for p in np.flatnonzero(seq == task.PROM):
        best = 0
        for e, t in enumerate(seq):
            if t not in task.STRENGTH:
                continue
            lo, hi = min(p, e), max(p, e)
            between = list(seq[lo + 1:hi])
            fwd = [m for m, s in enumerate(between) if s == task.FWD]
            rev = [m for m, s in enumerate(between) if s == task.REV]
            if not (fwd and rev and min(fwd) < max(rev)):
                best = max(best, task.STRENGTH[t])
        y[p] = best
    return y


class TestExtrusion(unittest.TestCase):
    def test_matches_brute_force(self):
        rng = np.random.default_rng(0)
        n = 9
        load, pl, pr = rng.uniform(0.05, 1, n), rng.uniform(0.05, 1, n), rng.uniform(0.05, 1, n)
        got = extrusion_log_contacts(
            torch.tensor(np.log(load)), torch.tensor(np.log(pl)), torch.tensor(np.log(pr))
        ).exp().numpy()
        np.fill_diagonal(got, 0)
        np.testing.assert_allclose(got, brute_extrusion(load, pl, pr), rtol=1e-9, atol=1e-12)

    def test_hard_barriers_reproduce_convergent_rule(self):
        """With '>' stopping left arms, '<' stopping right arms and uniform loading,
        contact is zero exactly when the task says the pair is insulated."""
        rng = np.random.default_rng(1)
        for _ in range(200):
            seq = task.sample_sequence(rng, 30)
            a = torch.tensor(np.where(seq == task.FWD, -50.0, 0.0))
            b = torch.tensor(np.where(seq == task.REV, -50.0, 0.0))
            d = extrusion_log_contacts(torch.zeros(30, dtype=torch.float64), a, b).exp().numpy()
            for i in range(30):
                for j in range(i + 1, 30):
                    between = list(seq[i + 1:j])
                    fwd = [m for m, s in enumerate(between) if s == task.FWD]
                    rev = [m for m, s in enumerate(between) if s == task.REV]
                    insulated = bool(fwd and rev and min(fwd) < max(rev))
                    self.assertEqual(d[i, j] < 1e-15, insulated, (seq, i, j))

    def test_one_arm_matches_definition(self):
        rng = np.random.default_rng(2)
        n = 7
        pl, pr = rng.uniform(0.1, 1, n), rng.uniform(0.1, 1, n)
        got = insulation_log_contacts(torch.tensor(np.log(pl)), torch.tensor(np.log(pr))).exp().numpy()
        for i in range(n):
            for j in range(n):
                if j > i:
                    want = np.prod(pr[i + 1:j])
                elif j < i:
                    want = np.prod(pl[j + 1:i])
                else:
                    want = 0.0
                self.assertAlmostEqual(got[i, j], want, places=12)


class TestTask(unittest.TestCase):
    def test_labels_match_brute_force(self):
        rng = np.random.default_rng(3)
        for _ in range(300):
            seq = task.sample_sequence(rng, int(rng.integers(5, 40)))
            np.testing.assert_array_equal(task.label_sequence(seq), brute_label(seq))


class TestModels(unittest.TestCase):
    def test_all_mixers_run_long_without_nan(self):
        tokens = torch.randint(0, 7, (2, 512))
        for name in MIXERS:
            model = SequenceModel(name)
            out = model(tokens)
            self.assertEqual(out.shape, (2, 512, 4))
            out.square().mean().backward()
            for p in model.parameters():
                if p.grad is not None:
                    self.assertTrue(torch.isfinite(p.grad).all(), name)


if __name__ == "__main__":
    unittest.main()
