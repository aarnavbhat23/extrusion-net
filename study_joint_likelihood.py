"""Frozen reanalysis of six existing occupancy experiments, retaining failures."""
import hashlib
import json
from pathlib import Path
import numpy as np
from joint_likelihood import split_counts, fit_training, log_likelihood, search_witnesses

ROOT = Path(__file__).resolve().parent


def main():
    source_path = ROOT / 'results/continuous_occupancy_study.json'
    old_path = ROOT / 'results/occupancy_ambiguity_diagnostic.json'
    source = json.loads(source_path.read_text())
    old = json.loads(old_path.read_text())
    rng = np.random.default_rng(20261028)
    rows = []
    for data, previous in zip(source['rows'], old['rows'], strict=True):
        assert data['parameters'] == previous['parameters'] and data['repeat'] == previous['repeat']
        splits, training, heldout = [], [], []
        for observation in data['observations']:
            train, test = split_counts(observation['counts'], rng)
            gamma = observation['gamma']
            splits.append({'gamma': gamma, 'training': train.tolist(), 'heldout': test.tolist()})
            training.append((gamma, train)); heldout.append((gamma, test))
        fitted = fit_training(training)
        threshold = log_likelihood(fitted['parameters'], heldout) + np.log(.05)
        witnesses = search_witnesses(heldout, fitted['parameters'], threshold)
        old_points = []
        for key in ('lower_absence_witness', 'upper_absence_witness'):
            p = previous[key]['parameters']
            slack = log_likelihood(p, heldout) - threshold
            old_points.append({'parameters': p, 'slack': float(slack), 'retained': bool(slack >= 0)})
        row = {'parameters': data['parameters'], 'repeat': data['repeat'], 'splits': splits,
               'fit': fitted, 'threshold': float(threshold), 'old_witnesses': old_points,
               'generating_slack': log_likelihood(data['parameters'], heldout) - threshold,
               'new_witnesses': witnesses}
        rows.append(row)
        print('case', len(rows), 'old retained', sum(p['retained'] for p in old_points),
              'new span', witnesses['attained_span'], 'generating slack', row['generating_slack'], flush=True)
    report = {'protocol_commit': '36acd46', 'seed': 20261028, 'alpha': .05,
              'scope': 'Development reanalysis; local inner witnesses, no global width certificate or biological validation.',
              'source_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'old_witness_sha256': hashlib.sha256(old_path.read_bytes()).hexdigest(), 'rows': rows}
    (ROOT / 'results/joint_likelihood.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
