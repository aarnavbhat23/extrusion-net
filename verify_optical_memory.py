"""Replay frozen count outcomes without importing the study's estimator."""
import json
from pathlib import Path
import numpy as np


def verify(root):
    result = json.loads((root/'results/optical_memory.json').read_text())
    checked = 0
    with np.load(root/'results/optical_memory_counts.npz', allow_pickle=False) as saved:
        assert len(saved.files) == len(result['settings'])
        for row in result['settings']:
            a, b, n, name = row['a'], row['b'], row['n_per_arm'], row['method']
            key = f'a{a}_b{b}_n{n}_{name}'
            counts = saved[key]
            assert counts.shape == (2000, 2) and np.all((counts >= 0)&(counts <= n))
            contrast = counts[:, 0]/n-counts[:, 1]/n
            radius = np.sqrt(np.log(40)/n)
            if name.startswith('primary_'):
                rejected = abs(contrast) > radius
            else:
                # Extremum of division by a positive interval depends on sign.
                left, right = contrast-radius-.02, contrast+radius+.02
                lower = np.maximum(-1., np.where(left < 0, left/.6, left/.8))
                upper = np.minimum(1., np.where(right < 0, right/.8, right/.6))
                rejected = lower > 0
                covered = (lower <= row['true_target'])&(row['true_target'] <= upper)
                assert abs(covered.mean()-row['coverage']) < 1e-12
            assert abs(rejected.mean()-row['rejection_rate']) < 1e-12
            checked += len(counts)
    verification = dict(status='passed', experiment_replicates_replayed=checked,
                        inference_refit=False, independent_biological_validation=False)
    (root/'results/optical_memory_verification.json').write_text(json.dumps(verification, indent=2)+'\n')
    print(json.dumps(verification, indent=2))


if __name__ == '__main__':
    verify(Path(__file__).resolve().parent)
