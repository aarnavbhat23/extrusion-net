"""Independent high-precision check of the explicit one-clock inverse."""
from decimal import Decimal, localcontext
import json
from pathlib import Path
from verify_occupancy_witnesses import decimal_snapshot

ROOT = Path(__file__).resolve().parent


def recover(probabilities, gamma, zero_tolerance=Decimal('1e-55')):
    """Algebraic audit on high-precision model probabilities, not noisy estimates.

    Values below zero_tolerance are unresolved numerically; the theorem has no
    such threshold. The audit cases keep positive capture well above it.
    """
    p = probabilities; g = Decimal.from_float(float(gamma))
    r = g/p[0]; mu = r-2-g
    pl = 1-r*r*p[4*9]/p[2*9]; pr = 1-r*r*p[4]/p[2]
    if max(abs(pl), abs(pr)) < zero_tolerance:
        return [mu, pl, pr, None]
    capture, near, before = (pl, p[3*9], p[2*9]) if pl >= pr else (pr, p[3], p[2])
    residual = near/before-(1-capture)/r
    h = (capture/residual-1-g)/mu
    return [mu, pl, pr, h]


def audit(parameters, gamma):
    probabilities = decimal_snapshot(parameters, gamma)
    estimated = recover(probabilities, gamma)
    errors = [abs(a-Decimal.from_float(float(b))) for a, b in zip(estimated, parameters) if a is not None]
    assert max(errors) < Decimal('1e-50')
    no_capture = parameters[1] == parameters[2] == 0
    assert (estimated[-1] is None) == no_capture
    if no_capture:
        alternative = list(parameters); alternative[-1] = 2.
        assert max(abs(a-b) for a, b in zip(decimal_snapshot(parameters, gamma, 2.), decimal_snapshot(alternative, gamma, 2.))) < Decimal('1e-55')
    return {'parameters': list(parameters), 'gamma': gamma,
            'recovered': [str(v) if v is not None else None for v in estimated],
            'maximum_error': str(max(errors)), 'h_irrelevant': no_capture}


def main():
    data = json.loads((ROOT/'results/continuous_occupancy_study.json').read_text())
    # Boundary cases supplement the existing mechanisms, with no data fitting.
    boundaries = ((.03, 0., 0., .15), (.3, 0., .7, 1.2), (.7, 1., 0., .5),
                  (1.5, 1., 1., 1.2), (.4, .000001, 0., .7))
    checked, diagnostics = [], []
    with localcontext() as context:
        context.prec = 70
        for p in boundaries:
            for gamma in (.1, .3):
                checked.append(audit(p, gamma))
        for row in data['rows']:
            for observation in row['observations']:
                gamma = observation['gamma']; p = row['parameters']
                checked.append(audit(p, gamma))
                probabilities = decimal_snapshot(p, gamma)
                n = sum(observation['counts'])
                diagnostics.append({'parameters': p, 'repeat': row['repeat'], 'gamma': gamma,
                                    'expected_beyond_left_count': float(n*probabilities[36]),
                                    'expected_beyond_right_count': float(n*probabilities[4]),
                                    'observed_beyond_left_count': observation['counts'][36],
                                    'observed_beyond_right_count': observation['counts'][4]})
    (ROOT/'results/capture_identifiability.json').write_text(json.dumps({
        'decimal_digits': 70, 'scope': 'Algebraic model audit and post-screen conditioning diagnosis; not a noisy estimator.',
        'checks': checked, 'existing_data_diagnostics': diagnostics}, indent=2)+'\n')
    print('Verified', len(checked), 'one-clock inverses including irrelevant-parameter boundary.')


if __name__ == '__main__':
    main()
