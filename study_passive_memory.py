"""Algebraic falsification of memory-attribution ideas; no biological fit."""
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm


def response(rates, initial, times):
    rates = np.asarray(rates, dtype=float)
    return np.exp(-np.outer(times, rates)) @ np.asarray(initial, dtype=float)


def passive_envelope(times, minimum_rate, initial_norm_bound, observation_norm):
    return 2*observation_norm*initial_norm_bound*np.exp(-minimum_rate*np.asarray(times))


def screen():
    rows = []
    for slow, fast in itertools.product([.1, .3, 1.], [1.5, 3., 10.]):
        peak_time = np.log(fast/slow)/(fast-slow)
        times = np.array([0., peak_time, 2*peak_time])
        means = response([slow, fast], [1., -1.], times)
        a, b = np.exp(-np.array([slow, fast])*peak_time)
        recurrence = np.array([a*b, -(a+b), 1.])
        offset_means = means+2.
        extra_mode = means+.4*np.exp(-.37*times)
        temporal_matrix = np.array([[1., 1.], [a, b]])
        rows.append(dict(slow=slow, fast=fast, times=times.tolist(),
                         preparation_A_means=means.tolist(), preparation_B_means=[0., 0., 0.],
                         both_initial_visible_variances=2., both_later_visible_variances=2.,
                         passive_mean_difference_at_peak=float(means[1]),
                         recurrence_residual=float(recurrence @ means),
                         constant_offset_residual=float(recurrence @ offset_means),
                         extra_passive_mode_residual=float(recurrence @ extra_mode),
                         temporal_observability_determinant=float(np.linalg.det(temporal_matrix)),
                         two_time_initial_mean_recovery=np.linalg.solve(temporal_matrix, means[:2]).tolist(),
                         independent_matrix_exponential_means=[float(np.ones(2) @ expm(-np.diag([slow, fast])*t) @ np.array([1., -1.])) for t in times],
                         passive_envelope=passive_envelope(times, slow, np.sqrt(2), np.sqrt(2)).tolist()))
    result = dict(scope='two independent passive Ornstein-Uhlenbeck modes, signed one-dimensional observation',
                  initial_mode_covariance='identity in both preparations; stationary covariance under noise sqrt(2*A)',
                  cases=rows, no_post_preparation_active_force=True, biological_validation=False,
                  new_architecture=False, conclusion='same present observation law does not imply same passive future law')
    destination = Path(__file__).resolve().parent/'results/passive_memory.json'
    destination.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(cases=len(rows), passive_peak_differences=[r['passive_mean_difference_at_peak'] for r in rows],
                          maximum_null_recurrence_residual=max(abs(r['recurrence_residual']) for r in rows),
                          minimum_offset_artifact=min(abs(r['constant_offset_residual']) for r in rows),
                          minimum_extra_mode_artifact=min(abs(r['extra_passive_mode_residual']) for r in rows)), indent=2))


if __name__ == '__main__':
    screen()
