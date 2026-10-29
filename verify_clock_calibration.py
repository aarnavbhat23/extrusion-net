"""Replay the saved study with a separate sparse-generator implementation.

Does not import the fitter, data generator, flow recurrence or interval code.
This checks computational consistency, not model adequacy or optimizer optimality.
"""
from functools import lru_cache
import json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parent
MECHANISMS = ((1.15,.85,.27,.35,.65,.45), (.75,1.25,.5,.6,.3,.9), (1.3,.7,.7,.45,.55,.6))
STRATEGIES = {'two_clocks': (.1,.3), 'close_clocks': (.28,.3), 'full_retention': (.3,)}


@lru_cache(None)
def snapshot(parameters, gamma, factor=1.):
    vl, vr, mu, pl, pr, h = parameters
    row, col, data = [], [], []
    def index(a,b,cl,cr): return 4*(9*a+b)+2*cl+cr
    for a in range(9):
        for b in range(9):
            for cl in range(2):
                for cr in range(2):
                    source = index(a,b,cl,cr)
                    outgoing = []
                    if a < 8 and cl == 0:
                        capture = pl if a == 2 else 0.
                        outgoing.extend(((index(a+1,b,0,cr),vl*(1-capture)),
                                         (index(a+1,b,1,cr),vl*capture)))
                    if b < 8 and cr == 0:
                        capture = pr if b == 2 else 0.
                        outgoing.extend(((index(a,b+1,cl,0),vr*(1-capture)),
                                         (index(a,b+1,cl,1),vr*capture)))
                    row.append(source); col.append(source)
                    data.append(gamma+factor*mu*(h if cl or cr else 1.)+sum(v for _,v in outgoing))
                    for target, rate in outgoing:
                        row.append(target); col.append(source); data.append(-rate)
    operator = coo_matrix((data,(row,col)),shape=(324,324)).tocsc()
    rhs = np.zeros(324); rhs[0] = gamma
    latent = spsolve(operator,rhs)
    present = latent.reshape(81,4).sum(axis=1)
    return np.r_[present,1-present.sum()]


def speed_interval(observations):
    intervals = []
    for obs in observations:
        counts = obs['counts']; successes = counts[0]; trials = successes+counts[9]
        if trials:
            ci = binomtest(successes,trials).proportion_ci(confidence_level=.975,method='exact')
            lower, upper = ci.low, ci.high
        else:
            lower, upper = 0.,1.
        intervals.append((lower/(1-lower),upper/(1-upper) if upper < 1 else np.inf))
    lower = intervals[1][0]-intervals[0][1]
    upper = intervals[1][1]-intervals[0][0]
    if upper <= 0: return None
    gap = observations[1]['gamma']-observations[0]['gamma']
    return [gap/upper, gap/lower if lower > 0 else np.inf]


def main():
    report = json.loads((ROOT/'results/clock_calibration.json').read_text())
    assert report['complete'] and report['protocol_commit'] == '60cf44c'
    assert report['seed'] == 20261029 and len(report['rows']) == 72
    rng = np.random.default_rng(20261029); cursor = 0
    for case, parameters in enumerate(MECHANISMS):
        target = snapshot(parameters,.3,2.)
        for count in (5000,20000):
            for repeat in range(4):
                for method, gammas in STRATEGIES.items():
                    r = report['rows'][cursor]; cursor += 1
                    assert (r['case'],r['initiated_episodes'],r['repeat'],r['method']) == (case,count,repeat,method)
                    assert tuple(r['parameters']) == parameters
                    conditional = method != 'full_retention'
                    assert r['conditional_likelihood'] == conditional
                    assert len(r['observations']) == len(gammas)
                    for obs, gamma in zip(r['observations'],gammas):
                        assert obs['gamma'] == gamma
                        generated = rng.multinomial(count//len(gammas),snapshot(parameters,gamma))
                        np.testing.assert_array_equal(obs['counts'],generated)
                    assert r['retained_episodes'] == sum(sum(o['counts'][:-1]) for o in r['observations'])
                    np.testing.assert_allclose(r['expected_sum_clock_duration'],sum(count/len(gammas)/g for g in gammas))
                    assert len(r['fit']['starts']) == 3
                    for fitted in r['fit']['starts']+[r['fit']]:
                        likelihood = 0.
                        pfit = np.array(fitted['parameters'])
                        assert np.all(pfit >= [.2,.2,.03,.001,.001,.1])
                        assert np.all(pfit <= [3,3,2,.999,.999,1.5])
                        for obs in r['observations']:
                            p = snapshot(tuple(pfit),obs['gamma'])
                            counts = np.array(obs['counts'])
                            if conditional:
                                p = p[:-1]/p[:-1].sum(); counts = counts[:-1]
                            mask = counts > 0
                            likelihood += counts[mask] @ np.log(p[mask])
                        np.testing.assert_allclose(likelihood,fitted['log_likelihood'],atol=1e-7,rtol=0)
                    assert r['fit']['log_likelihood'] == max(s['log_likelihood'] for s in r['fit']['starts'])
                    prediction = snapshot(tuple(r['fit']['parameters']),.3,2.)
                    np.testing.assert_allclose(r['prediction'],prediction,atol=1e-12,rtol=0)
                    np.testing.assert_allclose(r['target'],target,atol=1e-12,rtol=0)
                    np.testing.assert_allclose(r['absence_absolute_error'],abs(prediction[-1]-target[-1]),atol=1e-12,rtol=0)
                    np.testing.assert_allclose(r['target_total_variation'],np.abs(prediction-target).sum()/2,atol=1e-12,rtol=0)
                    if conditional:
                        ci = speed_interval(r['observations']); saved = r['speed_interval']
                        assert saved['empty'] == (ci is None)
                        if ci is not None:
                            encoded = [np.inf if v is None else v for v in saved['bounds']]
                            np.testing.assert_allclose(encoded,ci,atol=1e-8,rtol=1e-8)
                        assert saved['covers_generating_speed'] == (ci is not None and ci[0] <= parameters[0] <= ci[1])
                    else:
                        assert r['speed_interval'] is None
    rows = report['rows']; summary = []
    for n in (5000,20000):
        for method in STRATEGIES:
            selected = [r for r in rows if r['initiated_episodes'] == n and r['method'] == method]
            medians = [float(np.median([r['absence_absolute_error'] for r in selected if r['case'] == c])) for c in range(3)]
            summary.append({'initiated_episodes': n,'method': method,
                            'median_absence_error': float(np.median([r['absence_absolute_error'] for r in selected])),
                            'median_total_variation': float(np.median([r['target_total_variation'] for r in selected])),
                            'per_mechanism_median_absence_error': medians,
                            'successful_selected_fits': sum(r['fit']['success'] for r in selected),
                            'unbounded_speed_intervals': sum(r['speed_interval'] is not None and r['speed_interval']['bounds'] is not None and r['speed_interval']['bounds'][1] is None for r in selected)})
    primary = next(s for s in summary if s['initiated_episodes'] == 20000 and s['method'] == 'two_clocks')
    passed = sum(v <= .02 for v in primary['per_mechanism_median_absence_error']) >= 2
    output = {'verified_rows': cursor,'independent_sparse_replay': True,'summary': summary,'continuation_gate_passed': passed,
              'limitations': 'No refit or global optimum certification; shared assumed family and floating-point libraries; no biological validation.'}
    (ROOT/'results/clock_calibration_verification.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps(output,indent=2))


if __name__ == '__main__': main()
