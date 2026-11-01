"""Frozen loading-drift comparison, keeping common-loading failures visible."""
import gzip
import itertools
import json
from pathlib import Path
import numpy as np
from loading_target import binomial_interval, drift_bounds
from study_loading_target import wilson

ROOT = Path(__file__).resolve().parent
MECHANISMS = [(.2, .2, .75), (.4, .1, .85), (.65, .12, .95)]


def run():
    rng = np.random.default_rng(2026110103)
    rows, draws, certificates = [], [], []
    for index, ratio in itertools.product(range(3), [1., .9, .8]):
        mu, delta, ell = MECHANISMS[index]
        means = [ell*.3/(.3+mu+delta), ratio*ell*.3/(.3+2*mu+delta)]
        target = 1-ell*.3/(.3+3*mu+delta)
        contracts = dict(common_assumption=(1., 1.), broad_ratio=(.8, 1.),
                         calibrated_ratio=(max(.8, ratio-.02), min(1., ratio+.02)), exact_ratio=(ratio, ratio))
        for repeat in range(100):
            x = rng.binomial(12000, means)
            key = f'm{index}_k{ratio}_r{repeat}'
            draws.append(dict(key=key, mechanism=index, ratio=ratio, counts=x.tolist()))
            intervals = [binomial_interval(int(k), 12000, .025) for k in x]
            for method, contract in contracts.items():
                result = drift_bounds(intervals, ratio=contract)
                exact = drift_bounds([(p, p) for p in means], ratio=contract)
                low, high = result['interval']
                rows.append(dict(key=key, mechanism=index, ratio=ratio, method=method, target=target,
                                 interval=[float(low), float(high)], covered=bool(low <= target <= high),
                                 width=float(high-low), exact_identified_width=None if exact['status'] != 'ok' else float(exact['interval'][1]-exact['interval'][0]),
                                 exact_status=exact['status'], status=result['status']))
                certificates.append(dict(key=key, method=method, ratio_contract=contract, **result))
    summaries = []
    for index, ratio, method in itertools.product(range(3), [1., .9, .8], contracts):
        selected = [r for r in rows if (r['mechanism'], r['ratio'], r['method']) == (index, ratio, method)]
        covered = sum(r['covered'] for r in selected)
        summaries.append(dict(mechanism=index, ratio=ratio, method=method, coverage=covered/100,
                              coverage_wilson=wilson(covered, 100), median_width=float(np.median([r['width'] for r in selected])),
                              exact_identified_width=selected[0]['exact_identified_width'], exact_status=selected[0]['exact_status'],
                              incompatible=sum(r['status'] != 'ok' for r in selected)))
    gates = {}
    for method in contracts:
        selected = [r for r in summaries if r['method'] == method]
        gates[method] = dict(coverage=all(r['coverage_wilson'][0] >= .90 for r in selected),
                             precision=all(sum(r['median_width'] <= .10 for r in selected if r['ratio'] == k) >= 2 for k in [1., .9, .8]))
    result = dict(protocol_commit='868f4da', seed=2026110103, mechanisms=MECHANISMS,
                  settings=summaries, gates=gates, experiments=rows, draws=draws,
                  accurate_loading_reference_is_extra_information=True, loading_reference_measured=False,
                  biological_validation=False, unchanged_target_loading_required=True)
    (ROOT/'results/loading_drift.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    (ROOT/'results/loading_drift_certificates.json.gz').write_bytes(gzip.compress(json.dumps(certificates, separators=(',', ':'), allow_nan=False).encode(), mtime=0))
    print(json.dumps(dict(gates=gates, summary=summaries), indent=2))


if __name__ == '__main__':
    run()
