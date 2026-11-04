"""Frozen finite-data comparison with complete outcomes and resource accounting."""
import json
from pathlib import Path
import time
import numpy as np
from clock_calibration import generator_snapshot, flow_snapshot, fit, speed_interval

ROOT = Path(__file__).resolve().parent
MECHANISMS = ((1.15,.85,.27,.35,.65,.45), (.75,1.25,.5,.6,.3,.9), (1.3,.7,.7,.45,.55,.6))
STRATEGIES = {'two_clocks': (.1, .3), 'close_clocks': (.28, .3), 'full_retention': (.3,)}


def main():
    rng = np.random.default_rng(20261029); rows = []
    for case, parameters in enumerate(MECHANISMS):
        target = generator_snapshot(parameters, .3, 2.)
        for count in (5000, 20000):
            for repeat in range(4):
                for method, gammas in STRATEGIES.items():
                    conditional = method != 'full_retention'
                    observations, saved = [], []
                    for gamma in gammas:
                        probabilities = generator_snapshot(parameters, gamma)
                        np.testing.assert_allclose(probabilities, flow_snapshot(parameters, gamma), atol=1e-13, rtol=0)
                        counts = rng.multinomial(count//len(gammas), probabilities)
                        observations.append((gamma, counts[:-1] if conditional else counts))
                        saved.append({'gamma': gamma, 'counts': counts.tolist()})
                    started = time.perf_counter(); fitted = fit(observations, conditional)
                    seconds = time.perf_counter()-started
                    prediction = flow_snapshot(fitted['parameters'], .3, 2.)
                    interval = speed_interval(observations) if conditional else None
                    if interval is not None:
                        bounds = interval['bounds']
                        interval['covers_generating_speed'] = bool(bounds is not None and bounds[0] <= parameters[0]
                                                                  and (bounds[1] is None or parameters[0] <= bounds[1]))
                    row = {'case': case, 'parameters': parameters, 'initiated_episodes': count, 'repeat': repeat,
                           'method': method, 'observations': saved, 'conditional_likelihood': conditional,
                           'retained_episodes': sum(int(sum(x['counts'][:-1])) for x in saved),
                           'expected_sum_clock_duration': sum(count/len(gammas)/g for g in gammas),
                           'fit': fitted, 'fit_seconds': seconds, 'speed_interval': interval,
                           'prediction': prediction.tolist(), 'target': target.tolist(),
                           'target_total_variation': float(np.abs(prediction-target).sum()/2),
                           'absence_absolute_error': float(abs(prediction[-1]-target[-1]))}
                    rows.append(row)
                print('case', case, 'N', count, 'repeat', repeat,
                      [(x['method'], round(x['absence_absolute_error'], 5)) for x in rows[-3:]], flush=True)
                # Recoverable progress without treating partial output as completed evaluation.
                (ROOT/'results/clock_calibration_progress.json').write_text(json.dumps({'complete': False, 'rows': rows}, indent=2)+'\n')
    report = {'protocol_commit': '60cf44c', 'seed': 20261029, 'complete': True,
              'scope': 'Matched-family development simulations; no clinical validation, timing savings or prediction coverage guarantee.', 'rows': rows}
    (ROOT/'results/clock_calibration.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    (ROOT/'results/clock_calibration_progress.json').unlink()


if __name__ == '__main__':
    main()
