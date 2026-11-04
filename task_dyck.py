"""Generality test 1: nested brackets, nothing to do with genomes.

Tokens: '.' filler, '(' ')' '[' ']'. Sequences are prefixes of Dyck-2 words with
nesting depth at most MAX_DEPTH, so depth does not grow with length.

Label at every position: the type of the innermost bracket still open after
that token (0 none, 1 '(', 2 '['). After a closing bracket this requires
skipping every matched pair back to the right opener, a stack operation.
"""

import numpy as np

TOKENS = [".", "(", ")", "[", "]"]
DOT, OPEN_R, CLOSE_R, OPEN_S, CLOSE_S = range(5)
N_CLASSES = 3
IGNORE = -100
MAX_DEPTH = 6
CLOSER = {OPEN_R: CLOSE_R, OPEN_S: CLOSE_S}
LABEL = {OPEN_R: 1, OPEN_S: 2}


def sample_sequence(rng, n):
    stack, seq = [], []
    for _ in range(n):
        options, weights = [DOT], [0.4]
        if len(stack) < MAX_DEPTH:
            options += [OPEN_R, OPEN_S]
            weights += [0.15, 0.15]
        if stack:
            options.append(CLOSER[stack[-1]])
            weights.append(0.3)
        w = np.array(weights) / sum(weights)
        t = options[rng.choice(len(options), p=w)]
        if t in CLOSER:
            stack.append(t)
        elif t != DOT:
            stack.pop()
        seq.append(t)
    return np.array(seq)


def label_sequence(seq):
    stack, y = [], np.zeros(len(seq), dtype=np.int64)
    for i, t in enumerate(seq):
        if t in CLOSER:
            stack.append(t)
        elif t != DOT:
            stack.pop()
        y[i] = LABEL[stack[-1]] if stack else 0
    return y


def make_split(seed, n_seqs, lengths):
    rng = np.random.default_rng(seed)
    buckets = {}
    for _ in range(n_seqs):
        n = int(rng.choice(lengths))
        x = sample_sequence(rng, n)
        buckets.setdefault(n, ([], []))
        buckets[n][0].append(x)
        buckets[n][1].append(label_sequence(x))
    return {n: (np.stack(xs), np.stack(ys)) for n, (xs, ys) in buckets.items()}
