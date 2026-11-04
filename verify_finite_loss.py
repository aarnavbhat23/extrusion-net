"""Replay saved fits through the independently assembled full-state generator."""
import json
from pathlib import Path
import numpy as np
from finite_loss import reference, TIMES

ROOT = Path(__file__).resolve().parent


def log_likelihood(parameters, observations, method, partition):
    value = 0.
    for obs in observations:
        p = reference(tuple(parameters), obs['factor'], obs['gamma'])
        if obs['time'] is not None:
            p = p[TIMES.index(obs['time'])]
        counts = np.array(obs[partition + '_counts'])
        if method == 'three_normalized':
            p = p[:-1] / p[:-1].sum()
            counts = counts[:-1]
        mask = counts > 0
        value += counts[mask] @ np.log(p[mask])
    return float(value)


def main():
    data = json.loads((ROOT / 'results/finite_loss.json').read_text())
    rows = data['rows']
    assert data['complete'] and len(rows) == 90
    keys = {(r['stress'], r['case'], r['initiated_episodes'], r['repeat'], r['method']) for r in rows}
    assert len(keys) == 90
    rng = np.random.default_rng(data['seed'])
    max_likelihood_error = 0.
    for row in rows:
        observations, method = row['observations'], row['method']
        extra = .15 if row['stress'] == 'state_dropout' else 0.
        n = row['initiated_episodes'] // (2 * len(observations))
        for obs in observations:
            p = reference(tuple(row['parameters']), obs['actual_factor'], obs['gamma'], extra)
            if obs['time'] is not None:
                p = p[TIMES.index(obs['time'])]
            for partition in ('training', 'validation'):
                np.testing.assert_array_equal(rng.multinomial(n, p), obs[partition + '_counts'])
        assert sum(sum(o['training_counts']) + sum(o['validation_counts']) for o in observations) == row['initiated_episodes']
        assert len(row['fit']['starts']) == 3
        for fit in row['fit']['starts']:
            replay = log_likelihood(fit['parameters'], observations, method, 'training')
            error = abs(replay - fit['log_likelihood'])
            max_likelihood_error = max(max_likelihood_error, error)
            assert error < 1e-7
        assert row['fit']['log_likelihood'] == max(x['log_likelihood'] for x in row['fit']['starts'])
        for key, parameters in (('holdout_log_likelihood_at_fit', row['fit']['parameters']),
                                ('holdout_log_likelihood_at_base_parameter', row['parameters'])):
            replay = log_likelihood(parameters, observations, method, 'validation')
            max_likelihood_error = max(max_likelihood_error, abs(replay - row[key]))
            assert abs(replay - row[key]) < 1e-7
        np.testing.assert_allclose(row['holdout_threshold'], row['holdout_log_likelihood_at_fit'] + np.log(.05), atol=1e-9, rtol=0)
        assert row['base_parameter_in_split_set'] == (row['holdout_log_likelihood_at_base_parameter'] >= row['holdout_threshold'])
        target = reference(tuple(row['parameters']), 3., .3, extra)
        prediction = reference(tuple(row['fit']['parameters']), 3., .3)
        np.testing.assert_allclose(target, row['target'], atol=1e-12, rtol=0)
        np.testing.assert_allclose(prediction, row['prediction'], atol=1e-12, rtol=0)
        assert abs(row['absence_error'] - abs(prediction[-1] - target[-1])) < 1e-12
        assert abs(row['target_total_variation'] - abs(prediction - target).sum()/2) < 1e-12
    summary = []
    for stress in ('matched', 'turnover_error', 'state_dropout'):
        for n in ((6000, 24000) if stress == 'matched' else (24000,)):
            for method in ('three_normalized', 'absolute_retention', 'scheduled_timestamps'):
                selected = [r for r in rows if (r['stress'], r['initiated_episodes'], r['method']) == (stress, n, method)]
                medians = [float(np.median([r['absence_error'] for r in selected if r['case'] == case])) for case in range(3)]
                summary.append(dict(stress=stress, initiated_episodes=n, method=method, fits=len(selected),
                                    median_absence_error=float(np.median([r['absence_error'] for r in selected])),
                                    median_total_variation=float(np.median([r['target_total_variation'] for r in selected])),
                                    case_median_absence_errors=medians, passing_cases=sum(x <= .02 for x in medians),
                                    selected_optimizer_failures=sum(not r['fit']['success'] for r in selected),
                                    all_start_failures=sum(not s['success'] for r in selected for s in r['fit']['starts']),
                                    base_parameter_memberships=sum(r['base_parameter_in_split_set'] for r in selected),
                                    expected_sum_observation_duration=selected[0]['expected_sum_observation_duration']))
    report = dict(verified_rows=len(rows), verified_starts=270, maximum_likelihood_discrepancy=max_likelihood_error,
                  summary=summary, scope='Independent probability backend replay; same mathematical model, not biological validation.')
    (ROOT / 'results/finite_loss_verification.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
