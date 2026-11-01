"""Post-screen stronger timestamp control using certified loading slices."""
import gzip
import json
from pathlib import Path
import numpy as np
from loading_target import linear_extreme
from study_loading_target import wilson

ROOT = Path(__file__).resolve().parent


def refine(original, bins=16):
    if original['status'] != 'ok':
        return dict(interval=[0., 1.], status=original['status'], bins=[])
    problem = original['problem']
    lo = original['certificates'][0]['value']; hi = -original['certificates'][1]['value']
    if hi < lo and lo-hi < 1e-10:
        hi = lo
    edges = np.linspace(lo, hi, bins+1) if hi > lo+1e-12 else [lo, hi]
    lower, upper, saved = [], [], []
    for left, right in zip(edges[:-1], edges[1:]):
        bounds = list(problem['bounds']); bounds[0] = [float(left), float(right)]
        low = linear_extreme([0., 1., 0., 1.], problem['matrix'], problem['rhs'], bounds)
        high = linear_extreme([0., -1., 0., -1.], problem['matrix'], problem['rhs'], bounds)
        if low is None or high is None:
            continue
        lower.append(float(np.exp(left)/(1-high['value'])))
        upper.append(float(np.exp(right)/(1+low['value'])))
        saved.append(dict(left=float(left), right=float(right), problem=dict(matrix=problem['matrix'], rhs=problem['rhs'], bounds=bounds),
                          certificates=[dict(objective=[0., 1., 0., 1.], **low), dict(objective=[0., -1., 0., -1.], **high)]))
    if not saved:
        return dict(interval=[0., 1.], status='incompatible', bins=[])
    return dict(interval=[max(0., 1-max(upper)), min(1., 1-min(lower))], status='ok', bins=saved,
                one_sided_relaxation_slack_bound=float(1-np.exp(-(hi-lo)/bins)), sharp=False)


def run():
    source = json.loads((ROOT/'results/loading_target.json').read_text())
    originals = json.loads(gzip.decompress((ROOT/'results/loading_target_certificates.json.gz').read_bytes()))
    chosen = {r['key']:r for r in source['experiments'] if not r['stress'] and r['n'] == 24000 and r['method'] == 'timestamp_outer'}
    rows, certificates = [], []
    for original in originals:
        if original['method'] != 'timestamp_outer' or original['key'] not in chosen:
            continue
        old = chosen[original['key']]
        result = refine(original)
        low, high = result['interval']
        rows.append(dict(key=old['key'], mechanism=old['mechanism'], target=old['target'], interval=[low, high],
                         coverage=bool(low <= old['target'] <= high), width=high-low,
                         previous_width=old['width'], slack=result.get('one_sided_relaxation_slack_bound'), status=result['status']))
        certificates.append(dict(key=old['key'], **result))
    summaries = []
    for mechanism in range(3):
        subset = [r for r in rows if r['mechanism'] == mechanism]
        k = sum(r['coverage'] for r in subset)
        summaries.append(dict(mechanism=mechanism, coverage=k/100, coverage_wilson=wilson(k, 100),
                              median_width=float(np.median([r['width'] for r in subset])),
                              previous_median_width=float(np.median([r['previous_width'] for r in subset])),
                              maximum_one_sided_slack=max(r['slack'] for r in subset if r['slack'] is not None)))
    output = dict(protocol_commit='recorded_before_refinement', reused_development_counts=True, sharp=False,
                  settings=summaries, experiments=rows, unchanged_acquisition_times=True)
    (ROOT/'results/timestamp_refinement.json').write_text(json.dumps(output, indent=2)+'\n')
    (ROOT/'results/timestamp_refinement_certificates.json.gz').write_bytes(gzip.compress(json.dumps(certificates, separators=(',', ':')).encode(), mtime=0))
    print(json.dumps(summaries, indent=2))


if __name__ == '__main__':
    run()
