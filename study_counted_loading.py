"""Frozen assay-design screen with observed, uncertain loading references."""
import gzip
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from counted_loading import counted_bounds
from study_loading_drift import MECHANISMS
from study_loading_target import wilson

ROOT = Path(__file__).resolve().parent


def run():
    rng = np.random.default_rng(2026110105)
    records, draws, summaries = [], [], []
    for stress, index, ratio in itertools.product([False, True], range(3), [.8, .9, 1.]):
        mu, delta, ell = MECHANISMS[index]
        retention = [ell*.3/(.3+mu+delta), ratio*ell*.3/(.3+2*mu+delta)]
        target = 1-ell*.3/(.3+3*mu+delta)
        detection_ratio = .9 if stress else 1.
        for repeat in range(100):
            retained = rng.binomial(24000, retention).tolist()
            for budget in [500, 2000, 8000, 32000]:
                positive = rng.binomial(budget, [.9*ell, .9*detection_ratio*ratio*ell]).tolist()
                key = f's{stress}_m{index}_k{ratio}_r{repeat}_n{budget}'
                data = dict(key=key, stress=stress, mechanism=index, ratio=ratio,
                            budget=budget, repeat=repeat, retained=retained, reference_positive=positive)
                draws.append(data)
                contracts = {'common_sensitivity':(1., 1.)}
                if stress:
                    contracts['bounded_sensitivity'] = (.9, 1.)
                for method, contract in contracts.items():
                    result = counted_bounds(retained, [24000]*2, positive, [budget]*2, contract)
                    low, high = result['interval']
                    records.append(dict(**data, method=method, target=target, covered=low <= target <= high,
                                        width=high-low, **result))
    for stress, index, ratio, budget in itertools.product([False, True], range(3), [.8, .9, 1.], [500, 2000, 8000, 32000]):
        for method in (['common_sensitivity', 'bounded_sensitivity'] if stress else ['common_sensitivity']):
            rows = [r for r in records if (r['stress'], r['mechanism'], r['ratio'], r['budget'], r['method']) == (stress, index, ratio, budget, method)]
            covered = sum(r['covered'] for r in rows)
            summaries.append(dict(stress=stress, mechanism=index, ratio=ratio, budget=budget, method=method,
                                  covered=covered, repetitions=len(rows), coverage=covered/len(rows),
                                  coverage_wilson=wilson(covered, len(rows)),
                                  median_width=float(np.median([r['width'] for r in rows])),
                                  fallbacks=sum(r['status'] != 'ok' for r in rows)))
    largest = [r for r in summaries if not r['stress'] and r['budget'] == 32000]
    gates = dict(primary_coverage=all(r['coverage_wilson'][0] >= .90 for r in largest),
                 primary_width=all(r['median_width'] <= .10 for r in largest),
                 bounded_sensitivity_coverage=all(r['coverage_wilson'][0] >= .90 for r in summaries
                    if r['stress'] and r['method'] == 'bounded_sensitivity'))
    details = ROOT/'results/counted_loading_details.json.gz'
    details.write_bytes(gzip.compress(json.dumps(dict(draws=draws, records=records), separators=(',', ':'), allow_nan=False).encode(), mtime=0))
    result = dict(protocol_commit='ccca9f1', seed=2026110105, mechanisms=MECHANISMS,
                  retention_episodes_per_arm=24000, reference_episodes_per_arm=[500, 2000, 8000, 32000],
                  settings=summaries, gates=gates, simulated_reference_counts=True,
                  measured_biological_calibration=False, zero_false_positive_assumption=True,
                  details_file=details.name, details_sha256=hashlib.sha256(details.read_bytes()).hexdigest())
    (ROOT/'results/counted_loading.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(gates=gates, largest_reference_budget=largest), indent=2))


if __name__ == '__main__':
    run()
