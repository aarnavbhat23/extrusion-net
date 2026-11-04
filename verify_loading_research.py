"""Replay saved count intervals and certify LP extrema without an optimizer."""
import gzip
import hashlib
import json
from functools import lru_cache
from pathlib import Path
import numpy as np
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parent


@lru_cache(maxsize=None)
def exact_interval(k, n, error):
    ci = binomtest(k, n).proportion_ci(confidence_level=1-error, method='exact')
    return [ci.low, ci.high]


def certificate(problem, solution):
    a = np.array(problem['matrix']); b = np.array(problem['rhs'])
    x = np.array(solution['primal']); c = np.array(solution['objective'])
    dual = np.array(solution['inequality_dual'])
    lower = np.array(solution['lower_dual']); upper = np.array(solution['upper_dual'])
    assert max(a @ x-b) < 1e-7
    assert max(dual) <= 1e-7 and min(lower) >= -1e-7 and max(upper) <= 1e-7
    np.testing.assert_allclose(a.T @ dual+lower+upper, c, atol=1e-7)
    value = float(b @ dual)
    for i, (lo, hi) in enumerate(problem['bounds']):
        if lo is not None:
            assert x[i] >= lo-1e-7; value += lo*lower[i]
        else:
            assert abs(lower[i]) < 1e-7
        if hi is not None:
            assert x[i] <= hi+1e-7; value += hi*upper[i]
        else:
            assert abs(upper[i]) < 1e-7
    np.testing.assert_allclose([c @ x, value], [solution['value']]*2, atol=1e-7)


def snapshot_problem(intervals, ratio, factor):
    ul, uu = 1/intervals[0][1], 1/intervals[0][0]
    vl, vu = 1/intervals[1][1], 1/intervals[1][0]
    tl, tu = 2/(factor[1]-1), 2/(factor[0]-1)
    if ratio is None:
        return dict(matrix=[[1, -1, 0], [-tl, tl, -1], [tu, -tu, 1], [-2, 0, 1]],
                    rhs=[0, 0, 0, -2], bounds=[[max(ul, 1), uu], [vl, vu], [0, None]])
    return dict(matrix=[[0, ratio[0], -1, 0], [0, -ratio[1], 1, 0], [1, 0, -1, 0],
                        [-tl, 0, tl, -1], [tu, 0, -tu, 1], [-2, 0, 0, 1]],
                rhs=[0, 0, 0, 0, 0, -2], bounds=[[max(ul, 1), uu], [vl, vu], [0, None], [0, None]])


def check_problem(actual, expected):
    np.testing.assert_allclose(actual['matrix'], expected['matrix'], atol=1e-8)
    np.testing.assert_allclose(actual['rhs'], expected['rhs'], atol=1e-8)
    for got, want in zip(actual['bounds'], expected['bounds']):
        for g, w in zip(got, want):
            if w is None:
                assert g is None
            else:
                np.testing.assert_allclose(g, w, atol=1e-7)


def timestamp_problem(intervals):
    tl, tu = 2/1.2, 2/.8
    matrix = [[0, 0, tl, -1], [0, 0, -tu, 1], [0, -2, 0, 1]]; rhs = [0, 0, 0]
    for interval, row in zip(intervals, [[1, -1, 0, 0], [1, -2, 0, 0], [1, -1, -1, 0], [1, -2, -2, 0]]):
        matrix.append(row); rhs.append(np.log(interval[1]))
        if interval[0] > 0:
            matrix.append((-np.array(row)).tolist()); rhs.append(-np.log(interval[0]))
    return dict(matrix=matrix, rhs=rhs, bounds=[[np.log(.6), 0], [0, None], [0, None], [0, None]])


def verify_study(name, drift=False):
    study = json.loads((ROOT/f'results/{name}.json').read_text())
    payload = (ROOT/'results'/study['details_file']).read_bytes()
    assert hashlib.sha256(payload).hexdigest() == study['details_sha256']
    study.update(json.loads(gzip.decompress(payload)))
    saved = json.loads(gzip.decompress((ROOT/f'results/{name}_certificates.json.gz').read_bytes()))
    records = {(r['key'], r['method']):r for r in study['experiments']}
    draws = {r['key']:r for r in study['draws']}
    rng = np.random.default_rng(study['seed'])
    for data in study['draws']:
        mu, delta, ell = study['mechanisms'][data['mechanism']]
        ratio = data['ratio'] if drift else (.8 if data['stress'] else 1.)
        means = [ell*.3/(.3+mu+delta), ratio*ell*.3/(.3+2*mu+delta)]
        n = 24000 if drift else data['n']
        expected = rng.binomial(n//2, means)
        np.testing.assert_array_equal(expected, data['counts'] if drift else data['absolute_counts'])
        if not drift:
            time_means = [loading*np.exp(-rate*t) for loading, rate in [(ell, mu+delta), (ratio*ell, 2*mu+delta)]
                          for t in [1/.3, 2/.3]]
            np.testing.assert_array_equal(rng.binomial(n//4, time_means), data['timestamp_counts'])
    intervals = 0; lp_count = 0
    for result in saved:
        row = records[result['key'], result['method']]; data = draws[result['key']]
        if result['status'] != 'ok':
            assert row['interval'] == [0., 1.]
            intervals += 1; continue
        if drift:
            cp = [exact_interval(k, 12000, .025) for k in data['counts']]
            expected = snapshot_problem(cp, result['ratio_contract'], [1.8, 2.2])
        elif result['method'] == 'timestamp_outer':
            cp = [exact_interval(k, data['n']//4, .0125) for k in data['timestamp_counts']]
            expected = timestamp_problem(cp)
        else:
            cp = [exact_interval(k, data['n']//2, .025) for k in data['absolute_counts']]
            factor = [2, 2] if result['method'] == 'absolute_known_factor' else [1.8, 2.2]
            expected = snapshot_problem(cp, None, factor)
        check_problem(result['problem'], expected)
        for solution in result['certificates']:
            certificate(result['problem'], solution); lp_count += 1
        solutions = result['certificates']
        if not drift and result['method'] == 'timestamp_outer':
            lmin, lmax = solutions[0]['value'], -solutions[1]['value']
            smin, smax = solutions[2]['value'], -solutions[3]['value']
            estimate = [max(0., 1-np.exp(lmax)/(1+smin)), min(1., 1-np.exp(lmin)/(1+smax))]
        else:
            estimate = [max(0., 1-1/solutions[0]['value']), min(1., 1+1/solutions[1]['value'])]
            # Check physical endpoint mechanisms directly, not only LP values.
            for i, solution in enumerate(solutions):
                values = solution['primal']
                u, v = values[:2]
                j, w = (values[2], values[3]) if drift else (v, values[2])
                mu = .3*w/2; delta = .3*(u-1)-mu
                factor = 1+.3*(j-u)/mu if mu > 1e-12 else None
                ratio = j/v if drift else 1.
                assert mu >= -1e-8 and delta >= -1e-8
                if factor is not None:
                    assert 1.8-1e-7 <= factor <= 2.2+1e-7
                    means = [.3/(.3+mu+delta), ratio*.3/(.3+factor*mu+delta)]
                    for value, (lo, hi) in zip(means, cp):
                        assert lo-1e-7 <= value <= hi+1e-7
                    target = 1-.3/(.3+3*mu+delta)
                    np.testing.assert_allclose(target, row['interval'][i], atol=1e-7)
        np.testing.assert_allclose(estimate, row['interval'], atol=1e-7)
        assert row['covered'] == (row['interval'][0] <= row['target'] <= row['interval'][1])
        intervals += 1
    return dict(intervals=intervals, optimality_certificates=lp_count, regenerated_experiments=len(study['draws']))


def verify_refinement():
    study = json.loads((ROOT/'results/timestamp_refinement.json').read_text())
    saved = json.loads(gzip.decompress((ROOT/'results/timestamp_refinement_certificates.json.gz').read_bytes()))
    rows = {r['key']:r for r in study['experiments']}
    originals = json.loads(gzip.decompress((ROOT/'results/loading_target_certificates.json.gz').read_bytes()))
    originals = {r['key']:r for r in originals if r['method'] == 'timestamp_outer'}
    count = 0
    for result in saved:
        low, high = [], []
        original = originals[result['key']]
        np.testing.assert_allclose(result['bins'][0]['left'], original['certificates'][0]['value'], atol=1e-8)
        np.testing.assert_allclose(result['bins'][-1]['right'], -original['certificates'][1]['value'], atol=1e-8)
        for item in result['bins']:
            expected = dict(original['problem'])
            expected['bounds'] = list(expected['bounds'])
            expected['bounds'][0] = [item['left'], item['right']]
            check_problem(item['problem'], expected)
            for solution in item['certificates']:
                certificate(item['problem'], solution); count += 1
            low.append(np.exp(item['left'])/(1-item['certificates'][1]['value']))
            high.append(np.exp(item['right'])/(1+item['certificates'][0]['value']))
        interval = [max(0., 1-max(high)), min(1., 1-min(low))]
        np.testing.assert_allclose(interval, rows[result['key']]['interval'], atol=1e-8)
        assert all(abs(a['right']-b['left']) < 1e-10 for a, b in zip(result['bins'][:-1], result['bins'][1:]))
    return dict(intervals=len(saved), optimality_certificates=count, sharp=False)


def verify_full_generator():
    from finite_loss import reference
    count = 0
    for name in ['loading_target', 'loading_drift']:
        source = json.loads((ROOT/f'results/{name}.json').read_text())
        for mu, delta, ell in source['mechanisms']:
            parameters = (1.1, .9, mu, .4, .7, 1., delta)
            for factor in [1., 2., 3.]:
                full = ell*reference(parameters, factor, .3)[:-1].sum()
                np.testing.assert_allclose(full, ell*.3/(.3+factor*mu+delta), atol=1e-12)
                count += 1
    return count


if __name__ == '__main__':
    result = dict(verified=True, loading_target=verify_study('loading_target'), loading_drift=verify_study('loading_drift', True),
                  timestamp_refinement=verify_refinement(), full_generator_survival_checks=verify_full_generator(), independent_biological_validation=False)
    (ROOT/'results/loading_research_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
