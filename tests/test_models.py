import unittest

import numpy as np
import torch

import task
import task_dyck
import task_recall
from models import (
    MIXERS,
    SequenceModel,
    extrusion_linear_scan,
    extrusion_log_contacts,
    insulation_log_contacts,
)


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

    def test_linear_scan_matches_dense(self):
        """The O(n) recurrences give exactly sum_j D(i, j) v_j."""
        rng = np.random.default_rng(4)
        n, d = 23, 5
        load, pl, pr = rng.uniform(0.05, 1, n), rng.uniform(0.05, 1, n), rng.uniform(0.05, 1, n)
        v = rng.normal(size=(n, d))
        dense = extrusion_log_contacts(
            torch.tensor(np.log(load)), torch.tensor(np.log(pl)), torch.tensor(np.log(pr))
        ).exp().numpy()
        np.fill_diagonal(dense, 0)
        got = extrusion_linear_scan(*(torch.tensor(t) for t in (load, pl, pr, v))).numpy()
        np.testing.assert_allclose(got, dense @ v, rtol=1e-10, atol=1e-12)

    def test_off_diagonal_blocks_have_rank_two(self):
        """D is 2-quasiseparable: any block with all rows left of all columns has rank <= 2."""
        rng = np.random.default_rng(5)
        n = 40
        load, pl, pr = rng.uniform(0.05, 1, n), rng.uniform(0.3, 1, n), rng.uniform(0.3, 1, n)
        d = extrusion_log_contacts(
            torch.tensor(np.log(load)), torch.tensor(np.log(pl)), torch.tensor(np.log(pr))
        ).exp().numpy()
        for cut in (10, 20, 30):
            s = np.linalg.svd(d[:cut, cut:], compute_uv=False)
            self.assertLess(s[2] / s[0], 1e-10, cut)
            self.assertGreater(s[1] / s[0], 1e-6, cut)

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


class TestGeneralityTasks(unittest.TestCase):
    def test_dyck_labels_track_the_stack_and_respect_depth(self):
        rng = np.random.default_rng(6)
        for _ in range(200):
            seq = task_dyck.sample_sequence(rng, int(rng.integers(5, 80)))
            depth, stack = 0, []
            for t, y in zip(seq, task_dyck.label_sequence(seq)):
                if t in (task_dyck.OPEN_R, task_dyck.OPEN_S):
                    stack.append(t)
                elif t != task_dyck.DOT:
                    self.assertEqual(task_dyck.CLOSER[stack.pop()], t)  # always well matched
                self.assertLessEqual(len(stack), task_dyck.MAX_DEPTH)
                self.assertEqual(y, task_dyck.LABEL[stack[-1]] if stack else 0)

    def test_recall_queries_follow_their_key(self):
        rng = np.random.default_rng(7)
        for _ in range(200):
            seq, y = task_recall.sample(rng, int(rng.integers(24, 100)))
            pairs = {seq[i]: seq[i + 1] for i in range(len(seq) - 1)
                     if task_recall.KEY0 <= seq[i] < task_recall.QUERY0}
            self.assertEqual(len(pairs), sum(task_recall.KEY0 <= t < task_recall.QUERY0 for t in seq))
            for p in np.flatnonzero(y != task_recall.IGNORE):
                key = seq[p] - task_recall.QUERY0 + task_recall.KEY0
                self.assertEqual(pairs[key] - task_recall.VALUE0, y[p])
                self.assertLess(max(i for i, t in enumerate(seq) if t == key), p)


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
