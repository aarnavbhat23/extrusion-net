"""Frozen finite-data common-loading target screen with saved LP certificates."""
import gzip
import itertools
import json
from pathlib import Path
import numpy as np
from loading_target import binomial_interval, snapshot_bounds, timestamp_bounds
from scipy.stats import norm

ROOT = Path(__file__).resolve().parent
MECHANISMS = [(0.25, .1, .7), (.5, .15, .85), (.7, .05, .95)]


def wilson(k, n):
    z = norm.ppf(.975); p = k/n; den = 1+z*z/n
    center = (p+z*z/(2*n))/den
    radius = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [float(center-radius), float(center+radius)]


def run():
    rng = np.random.default_rng(2026110101)
    rows, certificates, draws = [], [], []
    for stress, n, index in itertools.product([False, True], [6000, 24000], range(3)):
        if stress and n != 24000:
            continue
        mu, delta, ell = MECHANISMS[index]
        treated_loading = ell*.8 if stress else ell
        means = [ell*.3/(.3+mu+delta), treated_loading*.3/(.3+2*mu+delta)]
        time_means = [loading*np.exp(-rate*t) for loading, rate in [(ell, mu+delta), (treated_loading, 2*mu+delta)]
                      for t in [1/.3, 2/.3]]
        target = 1-ell*.3/(.3+3*mu+delta)
        exact = snapshot_bounds([(p, p) for p in means])['interval']
        for repeat in range(100):
            x = rng.binomial(n//2, means)
            y = rng.binomial(n//4, time_means)
            key = f'stress{stress}_n{n}_m{index}_r{repeat}'
            draws.append(dict(key=key, stress=stress, n=n, mechanism=index, repeat=repeat,
                              absolute_counts=x.tolist(), timestamp_counts=y.tolist()))
            intervals = [binomial_interval(int(k), n//2, .025) for k in x]
            time_intervals = [binomial_interval(int(k), n//4, .0125) for k in y]
            methods = dict(absolute_uncertain=snapshot_bounds(intervals),
                           absolute_known_factor=snapshot_bounds(intervals, factor=(2., 2.)),
                           timestamp_outer=timestamp_bounds(time_intervals))
            r0, r1 = x/(n//2)
            assumed_a = .3/r0-.3 if r0 > 0 else np.inf
            assumed_c = .3/r1-.3 if r1 > 0 else np.inf
            naive = 1-.3/(.3+assumed_a+2*(assumed_c-assumed_a))
            for method, result in methods.items():
                low, high = result['interval']
                rows.append(dict(key=key, stress=stress, n=n, mechanism=index, repeat=repeat, method=method,
                                 target=target, interval=[float(low), float(high)], covered=bool(low <= target <= high),
                                 width=float(high-low), exact_identified_width=float(exact[1]-exact[0]), status=result['status'],
                                 fixed_loading_plugin_absolute_error=float(abs(naive-target))))
                certificates.append(dict(key=key, method=method, **result))
    summaries = []
    for stress, n, index, method in itertools.product([False, True], [6000, 24000], range(3),
                                                     ['absolute_uncertain', 'absolute_known_factor', 'timestamp_outer']):
        selected = [r for r in rows if (r['stress'], r['n'], r['mechanism'], r['method']) == (stress, n, index, method)]
        if not selected:
            continue
        covered = sum(r['covered'] for r in selected)
        summaries.append(dict(stress=stress, n=n, mechanism=index, method=method,
                              coverage=covered/len(selected), coverage_wilson=wilson(covered, len(selected)),
                              median_width=float(np.median([r['width'] for r in selected])),
                              exact_identified_width=selected[0]['exact_identified_width'],
                              incompatible=sum(r['status'] != 'ok' for r in selected),
                              median_fixed_loading_plugin_error=float(np.median([r['fixed_loading_plugin_absolute_error'] for r in selected]))))
    gates = {}
    for method in ['absolute_uncertain', 'absolute_known_factor', 'timestamp_outer']:
        chosen = [r for r in summaries if not r['stress'] and r['n'] == 24000 and r['method'] == method]
        gates[method] = dict(coverage=all(r['coverage_wilson'][0] >= .90 for r in chosen),
                             precision=sum(r['median_width'] <= .10 for r in chosen) >= 2)
    result = dict(protocol_commit='c2bc901', seed=2026110101, scope='h=1, shared unknown loading, binary survival target',
                  mechanisms=MECHANISMS, gamma=.3, actual_factor=2., allowed_factor=[1.8, 2.2], target_factor=3.,
                  settings=summaries, gates=gates, experiments=rows, draws=draws,
                  cost_ledger=dict(absolute_groups=2, timestamp_groups=4,
                                   absolute_expected_duration_per_episode=1/.3, timestamp_mean_duration_per_episode=1.5/.3,
                                   fixed_loading_plugin_is_bias_diagnostic=True, known_factor_is_extra_information=True),
                  biological_validation=False, general_hidden_generator_guarantee=False)
    (ROOT/'results/loading_target.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    (ROOT/'results/loading_target_certificates.json.gz').write_bytes(gzip.compress(json.dumps(certificates, separators=(',', ':'), allow_nan=False).encode(), mtime=0))
    print(json.dumps(dict(gates=gates, primary=[r for r in summaries if not r['stress'] and r['n'] == 24000]), indent=2))


if __name__ == '__main__':
    run()
