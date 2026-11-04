"""Independent 70-digit likelihood replay using the separate decimal forward flow."""
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import numpy as np
from verify_occupancy_witnesses import decimal_snapshot

ROOT = Path(__file__).resolve().parent


def likelihood(parameters, splits):
    total = Decimal(0)
    for split in splits:
        probabilities = decimal_snapshot(parameters, split['gamma'])
        for count, p in zip(split['heldout'], probabilities, strict=True):
            if count:
                if p == 0:
                    return Decimal('-Infinity')
                total += count * p.ln()
    return total


def main():
    report = json.loads((ROOT / 'results/joint_likelihood.json').read_text())
    source_path = ROOT / 'results/continuous_occupancy_study.json'
    old_path = ROOT / 'results/occupancy_ambiguity_diagnostic.json'
    assert hashlib.sha256(source_path.read_bytes()).hexdigest() == report['source_sha256']
    assert hashlib.sha256(old_path.read_bytes()).hexdigest() == report['old_witness_sha256']
    source = json.loads(source_path.read_text())
    rng = np.random.default_rng(report['seed'])
    verified = []
    with localcontext() as context:
        context.prec = 70
        for data, row in zip(source['rows'], report['rows'], strict=True):
            for observation, split in zip(data['observations'], row['splits'], strict=True):
                counts = np.array(observation['counts'])
                train = rng.multivariate_hypergeometric(counts, int(counts.sum())//2)
                assert train.tolist() == split['training']
                assert (counts-train).tolist() == split['heldout']
            threshold = likelihood(row['fit']['parameters'], row['splits']) + Decimal('.05').ln()
            assert abs(threshold-Decimal.from_float(row['threshold'])) < Decimal('1e-8')
            slacks, targets = [], []
            for name in ('lower', 'upper'):
                point = row['new_witnesses'][name]
                if point is None:
                    continue
                assert all(l <= v <= h for l, v, h in zip((.03, 0, 0, .15), point['parameters'], (1.5, 1, 1, 1.2)))
                slack = likelihood(point['parameters'], row['splits'])-threshold
                assert slack > Decimal('1e-6')
                assert abs(slack-Decimal.from_float(point['slack'])) < Decimal('1e-8')
                target = decimal_snapshot(point['parameters'], factor=2.)[-1]
                assert abs(target-Decimal.from_float(point['absence'])) < Decimal('1e-12')
                slacks.append(str(slack)); targets.append(target)
            old_retained = []
            for point in row['old_witnesses']:
                slack = likelihood(point['parameters'], row['splits'])-threshold
                assert abs(slack-Decimal.from_float(point['slack'])) < Decimal('1e-8')
                assert (slack >= 0) == point['retained']
                old_retained.append(slack >= 0)
            span = targets[1]-targets[0] if len(targets) == 2 else None
            if span is not None:
                assert abs(span-Decimal.from_float(row['new_witnesses']['attained_span'])) < Decimal('1e-12')
            generating_slack = likelihood(row['parameters'], row['splits'])-threshold
            assert abs(generating_slack-Decimal.from_float(row['generating_slack'])) < Decimal('1e-8')
            verified.append({'parameters': row['parameters'], 'repeat': row['repeat'],
                             'old_retained': old_retained, 'new_slacks': slacks,
                             'attained_span': str(span), 'exceeds_point02': span > Decimal('.02') if span is not None else None,
                             'generating_slack': str(generating_slack)})
    output = {'decimal_digits': 70, 'verified_splits': True, 'rows': verified}
    (ROOT / 'results/joint_likelihood_verification.json').write_text(json.dumps(output, indent=2) + '\n')
    print('Verified all splits, old/new witnesses and generating parameters at 70 digits.')


if __name__ == '__main__':
    main()
