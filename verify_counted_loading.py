"""Replay reference counts and check interval/LP certificates independently."""
import gzip
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from verify_loading_research import certificate, check_problem, exact_interval, snapshot_problem

ROOT = Path(__file__).resolve().parent


def verify():
    study = json.loads((ROOT/'results/counted_loading.json').read_text())
    payload = (ROOT/'results'/study['details_file']).read_bytes()
    assert hashlib.sha256(payload).hexdigest() == study['details_sha256']
    details = json.loads(gzip.decompress(payload))
    rng = np.random.default_rng(2026110105)
    generated = []
    for stress, index, ratio in itertools.product([False, True], range(3), [.8, .9, 1.]):
        mu, delta, ell = study['mechanisms'][index]
        for repeat in range(100):
            retained = rng.binomial(24000, [ell*.3/(.3+mu+delta), ratio*ell*.3/(.3+2*mu+delta)]).tolist()
            for budget in [500, 2000, 8000, 32000]:
                positive = rng.binomial(budget, [.9*ell, .9*(.9 if stress else 1)*ratio*ell]).tolist()
                generated.append(dict(key=f's{stress}_m{index}_k{ratio}_r{repeat}_n{budget}', stress=stress,
                                      mechanism=index, ratio=ratio, budget=budget, repeat=repeat,
                                      retained=retained, reference_positive=positive))
    assert generated == details['draws']
    count = 0
    for row in details['records']:
        retention = [exact_interval(k, 24000, .0125) for k in row['retained']]
        reference = [exact_interval(k, row['budget'], .0125) for k in row['reference_positive']]
        np.testing.assert_allclose(row['retention_intervals'], retention, atol=1e-9)
        np.testing.assert_allclose(row['reference_intervals'], reference, atol=1e-9)
        qlo, qhi = row['sensitivity_ratio']
        assert list(row['sensitivity_ratio']) == ([.9, 1.] if row['method'] == 'bounded_sensitivity' else [1., 1.])
        low, high = reference[1][0]/reference[0][1]/qhi, min(1., reference[1][1]/reference[0][0]/qlo)
        if low > high:
            assert row['status'] == 'incompatible_reference' and row['interval'] == [0., 1.]
            continue
        np.testing.assert_allclose(row['ratio_contract'], [low, high], atol=1e-8)
        if row['status'] == 'ok':
            check_problem(row['problem'], snapshot_problem(retention, [low, high], [1.8, 2.2]))
            for endpoint, solution in enumerate(row['certificates']):
                certificate(row['problem'], solution); count += 1
                u, v, j, w = solution['primal']
                mu, delta, kappa = .3*w/2, .3*(u-1)-.3*w/2, j/v
                assert mu >= -1e-7 and delta >= -1e-7
                assert low-1e-7 <= kappa <= high+1e-7
                if mu > 1e-10:
                    factor = 1+.3*(j-u)/mu
                    assert 1.8-1e-7 <= factor <= 2.2+1e-7
                    for p, (lo, hi) in zip([.3/(.3+mu+delta), kappa*.3/(.3+factor*mu+delta)], retention):
                        assert lo-1e-7 <= p <= hi+1e-7
                target = 1-.3/(.3+3*mu+delta)
                np.testing.assert_allclose(target, row['interval'][endpoint], atol=1e-7)
        else:
            assert row['interval'] == [0., 1.]
        mu, delta, ell = study['mechanisms'][row['mechanism']]
        target = 1-ell*.3/(.3+3*mu+delta)
        np.testing.assert_allclose(row['target'], target, atol=1e-12)
        assert row['covered'] == (row['interval'][0] <= target <= row['interval'][1])
    for setting in study['settings']:
        rows = [r for r in details['records'] if all(r[k] == setting[k] for k in ['stress', 'mechanism', 'ratio', 'budget', 'method'])]
        assert len(rows) == setting['repetitions'] == 100
        assert sum(r['covered'] for r in rows) == setting['covered']
        np.testing.assert_allclose(np.median([r['width'] for r in rows]), setting['median_width'])
    output = dict(regenerated_reference_datasets=len(generated), interval_images=len(details['records']),
                  checked_optimality_certificates=count, optimizer_used_in_verification=False,
                  biological_validation=False)
    (ROOT/'results/counted_loading_verification.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    verify()
