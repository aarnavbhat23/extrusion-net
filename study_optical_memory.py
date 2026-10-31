"""Frozen biological/optical attribution screen; synthetic, calibrated reporters."""
import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm


def transitions(a, b, duration=1.):
    return {'A': expm(np.array([[-a, a], [0., 0.]])*duration),
            'B': expm(np.array([[0., 0.], [b, -b]])*duration)}


def occupancy(word, a, b):
    state = np.array([1., 0.])
    matrices = transitions(a, b)
    for letter in word:
        state = state @ matrices[letter]
    return float(state[1])


def observation_probability(word, observations, a, b, emission):
    """Forward likelihood for a binary readout after each pulse."""
    state = np.array([1., 0.])
    matrices = transitions(a, b)
    emission = np.asarray(emission)
    for letter, observed in zip(word, observations, strict=True):
        state = state @ matrices[letter]
        state *= emission if observed else 1-emission
    return float(state.sum())


def contrast_interval(difference, radius=0., epsilon=.01, calibration=(.6, .8)):
    """Exact inversion of d=kappa*Delta+u at contrast level, |u|<=2 epsilon."""
    if not 0 < calibration[0] <= calibration[1] or epsilon < 0 or radius < 0:
        raise ValueError('Positive calibration and nonnegative uncertainty required')
    difference = np.asarray(difference)
    low = difference-radius-2*epsilon
    high = difference+radius+2*epsilon
    corners = np.stack([low/calibration[0], low/calibration[1],
                        high/calibration[0], high/calibration[1]])
    return np.maximum(-1., corners.min(axis=0)), np.minimum(1., corners.max(axis=0))


def wilson(successes, total):
    z = 1.959963984540054
    p = successes/total
    middle = (p+z*z/(2*total))/(1+z*z/total)
    radius = z*np.sqrt(p*(1-p)/total+z*z/(4*total*total))/(1+z*z/total)
    return [float(middle-radius), float(middle+radius)]


def law_audit():
    comparisons, maximum, sum_error = 0, 0., 0.
    for length in range(1, 5):
        for word in itertools.product('AB', repeat=length):
            probabilities = []
            for bits in itertools.product((0, 1), repeat=length):
                # Distinct biological interpretations, identical permitted primary law.
                biological = observation_probability(word, bits, .6, 1.2, [1., .4])
                optical = observation_probability(word, bits, .6, 1.2, [1., .4])
                maximum = max(maximum, abs(biological-optical))
                probabilities.append(biological)
                comparisons += 1
            sum_error = max(sum_error, abs(sum(probabilities)-1))
    return dict(observation_strings=comparisons, maximum_law_difference=maximum,
                maximum_probability_sum_error=sum_error)


def run(output):
    rng = np.random.default_rng(2026103001)
    rows, saved = [], {}
    repetitions, epsilon, kappa, f = 2000, .01, .7, .1
    for a, b, n in itertools.product([.2, .6, 1.5], [.3, 1.2], [500, 2000]):
        x_ab, x_ba = occupancy('AB', a, b), occupancy('BA', a, b)
        delta = x_ba-x_ab
        analytic = (1-np.exp(-a))*(1-np.exp(-b))
        if abs(delta-analytic) > 1e-12:
            raise AssertionError('Matrix and analytic occupancy disagree')
        radius = np.sqrt(np.log(2/.05)/n)
        # Probabilities are ordered so the expected contrast is positive.
        primary = (1-.6*x_ab, 1-.6*x_ba)
        alternatives = {
            'primary_biological': (*primary, delta, False),
            'primary_optical': (*primary, 0., False),
            'reporter_biological': (f+kappa*x_ba-epsilon, f+kappa*x_ab+epsilon, delta, True),
            'reporter_optical': (f+epsilon, f-epsilon, 0., True),
            'reporter_optical_crosstalk': (f+kappa*x_ba-epsilon, f+kappa*x_ab+epsilon, 0., True),
            'direct_A_biological': (f+kappa*occupancy('A', a, b)-epsilon, f+epsilon,
                                    occupancy('A', a, b), True),
            'direct_A_optical': (f+epsilon, f-epsilon, 0., True),
        }
        for name, (p_plus, p_minus, truth, calibrated) in alternatives.items():
            counts = rng.binomial(n, np.array([p_plus, p_minus]), size=(repetitions, 2))
            difference = (counts[:, 0]-counts[:, 1])/n
            if calibrated:
                lower, upper = contrast_interval(difference, radius)
                rejected = lower > 0
                covered = (lower <= truth)&(truth <= upper)
            else:
                rejected = abs(difference) > radius
                covered = None
            key = f'a{a}_b{b}_n{n}_{name}'
            saved[key] = counts.astype(np.int16)
            row = dict(a=a, b=b, n_per_arm=n, method=name, true_target=truth,
                       plus_probability=p_plus, minus_probability=p_minus,
                       rejection_rate=float(rejected.mean()),
                       rejection_wilson=wilson(int(rejected.sum()), repetitions),
                       coverage=None if covered is None else float(covered.mean()),
                       coverage_wilson=None if covered is None else wilson(int(covered.sum()), repetitions),
                       outside_calibration_contract=name.endswith('crosstalk'))
            rows.append(row)
    valid = [r for r in rows if r['method'].startswith(('reporter_', 'direct_'))
             and not r['outside_calibration_contract']]
    nulls = [r for r in valid if r['method'].endswith('_optical')]
    alternatives = [r for r in rows if r['method'] == 'reporter_biological' and r['n_per_arm'] == 2000]
    law = law_audit()
    gates = dict(law_equivalence=law['maximum_law_difference'] <= 1e-12,
                 minimum_coverage=min(r['coverage'] for r in valid) >= .94,
                 null_control=max(r['rejection_wilson'][1] for r in nulls) <= .075,
                 power_cases=sum(r['rejection_rate'] >= .80 for r in alternatives))
    result = dict(seed=2026103001, replicates_per_setting=repetitions,
                  protocol='docs/OPTICAL_MEMORY_PROTOCOL.md', protocol_commit='f006846',
                  synthetic=True, law_audit=law, settings=rows, gates=gates,
                  conditional_feasibility_pass=all([gates['law_equivalence'], gates['minimum_coverage'],
                                                    gates['null_control'], gates['power_cases'] >= 4]),
                  biological_validation=False, novel_method_established=False)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'optical_memory_counts.npz', **saved)
    (output/'optical_memory.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    with (output/'optical_memory_summary.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    print(json.dumps({k:v for k,v in result.items() if k != 'settings'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent/'results')
    run(parser.parse_args().output)
