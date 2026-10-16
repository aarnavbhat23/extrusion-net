"""Generality test 2: associative recall, which needs lookup by content.

Tokens: '.' filler, keys K0-K7, queries Q0-Q7, values V0-V7. Up to 8 distinct
key-value pairs (key immediately followed by its value) appear in the first 70%
of the sequence; queries appear in the last 30%. Label at each query Qk: the
value that followed key Kk. Every other position is ignored.

Extrusion has no query-key matching, so this is where it is expected to be
weakest. Attention is built for it.
"""

import numpy as np

N_KEYS = 8
TOKENS = ["."] + [f"K{i}" for i in range(N_KEYS)] + [f"Q{i}" for i in range(N_KEYS)] + [f"V{i}" for i in range(N_KEYS)]
DOT = 0
KEY0, QUERY0, VALUE0 = 1, 1 + N_KEYS, 1 + 2 * N_KEYS
N_CLASSES = N_KEYS
IGNORE = -100


def sample(rng, n):
    seq = np.full(n, DOT, dtype=np.int64)
    y = np.full(n, IGNORE, dtype=np.int64)
    split = int(0.7 * n)
    n_pairs = int(rng.integers(2, min(N_KEYS, split // 3) + 1))
    # Non-overlapping pair slots in the first part, at even offsets with gaps.
    slots = rng.choice(split // 2, size=n_pairs, replace=False) * 2
    keys = rng.choice(N_KEYS, size=n_pairs, replace=False)
    values = rng.integers(0, N_KEYS, size=n_pairs)
    for s, k, v in zip(slots, keys, values):
        seq[s], seq[s + 1] = KEY0 + k, VALUE0 + v
    n_queries = int(rng.integers(1, min(4, n - split) + 1))
    for p in rng.choice(np.arange(split, n), size=n_queries, replace=False):
        which = int(rng.integers(n_pairs))
        seq[p], y[p] = QUERY0 + keys[which], values[which]
    return seq, y


def make_split(seed, n_seqs, lengths):
    rng = np.random.default_rng(seed)
    buckets = {}
    for _ in range(n_seqs):
        n = int(rng.choice(lengths))
        x, y = sample(rng, n)
        buckets.setdefault(n, ([], []))
        buckets[n][0].append(x)
        buckets[n][1].append(y)
    return {n: (np.stack(xs), np.stack(ys)) for n, (xs, ys) in buckets.items()}
