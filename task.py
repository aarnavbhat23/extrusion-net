"""Toy genome task: which enhancer can switch on each promoter?

Tokens: '.' filler, E1/E2/E3 enhancers of strength 1-3, P promoter,
'>' forward CTCF site, '<' reverse CTCF site.

Rule (the convergent-CTCF rule of real chromatin): a promoter and an enhancer
are insulated when a '>' appears before a '<' strictly between them, because
that pair anchors a loop no extruding motor can cross. A lone site, or a
divergent '<' ... '>' pair, does not insulate.

Label at each promoter: strength of the strongest enhancer it can reach (0-3).
Every other position is ignored (-100).
"""

import numpy as np

TOKENS = [".", "E1", "E2", "E3", "P", ">", "<"]
DOT, E1, E2, E3, PROM, FWD, REV = range(7)
STRENGTH = {E1: 1, E2: 2, E3: 3}
# Tuned so no label exceeds ~45% at any test length.
PROBS = np.array([0.63, 0.02, 0.02, 0.02, 0.07, 0.12, 0.12])
IGNORE = -100


def sample_sequence(rng, n):
    return rng.choice(len(TOKENS), size=n, p=PROBS)


def label_sequence(seq):
    n = len(seq)
    y = np.full(n, IGNORE, dtype=np.int64)
    for p in np.flatnonzero(seq == PROM):
        best = 0
        # Rightward: insulated once a '<' follows a '>' (both left of the enhancer).
        seen_fwd = blocked = False
        for m in range(p + 1, n):
            t = seq[m]
            if t in STRENGTH and not blocked:
                best = max(best, STRENGTH[t])
            if t == FWD:
                seen_fwd = True
            elif t == REV and seen_fwd:
                blocked = True
                break
        # Leftward: positions arrive in reverse order, so look for '>' after '<'.
        seen_rev = blocked = False
        for m in range(p - 1, -1, -1):
            t = seq[m]
            if t in STRENGTH and not blocked:
                best = max(best, STRENGTH[t])
            if t == REV:
                seen_rev = True
            elif t == FWD and seen_rev:
                blocked = True
                break
        y[p] = best
    return y


def make_split(seed, n_seqs, lengths):
    """Sequences grouped by length so batches never need padding."""
    rng = np.random.default_rng(seed)
    buckets = {}
    for _ in range(n_seqs):
        n = int(rng.choice(lengths))
        x = sample_sequence(rng, n)
        buckets.setdefault(n, ([], []))
        buckets[n][0].append(x)
        buckets[n][1].append(label_sequence(x))
    return {n: (np.stack(xs), np.stack(ys)) for n, (xs, ys) in buckets.items()}


def class_counts(split):
    counts = np.zeros(4, dtype=np.int64)
    for _, y in split.values():
        vals = y[y != IGNORE]
        counts += np.bincount(vals, minlength=4)
    return counts
