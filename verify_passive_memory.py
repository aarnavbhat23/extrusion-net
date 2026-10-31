"""Check saved passive counterexamples with matrix propagation."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm


def verify():
    root = Path(__file__).resolve().parent
    result = json.loads((root/'results/passive_memory.json').read_text())
    for row in result['cases']:
        a = np.diag([row['slow'], row['fast']])
        c, initial = np.ones(2), np.array([1., -1.])
        assert c @ initial == 0 and c @ np.eye(2) @ c == 2
        times = np.array(row['times'])
        means = np.array([c @ expm(-a*t) @ initial for t in times])
        np.testing.assert_allclose(means, row['preparation_A_means'], atol=1e-12)
        np.testing.assert_allclose(means, row['independent_matrix_exponential_means'], atol=1e-12)
        assert means[1] > 0 and np.all(abs(means) <= np.array(row['passive_envelope'])+1e-12)
        aa, bb = np.diag(expm(-a*times[1]))
        contrast = np.array([aa*bb, -(aa+bb), 1.])
        np.testing.assert_allclose(contrast @ means, row['recurrence_residual'], atol=1e-12)
        np.testing.assert_allclose(2*contrast.sum(), row['constant_offset_residual'], atol=1e-12)
        third = .4*np.array([expm(np.array([[-.37*t]]))[0, 0] for t in times])
        np.testing.assert_allclose(contrast @ third, row['extra_passive_mode_residual'], atol=1e-12)
        np.testing.assert_allclose(row['two_time_initial_mean_recovery'], initial, atol=1e-12)
        assert abs(row['constant_offset_residual']) > .01
        assert abs(row['extra_passive_mode_residual']) > .001
    checked = dict(verified=True, passive_cases=len(result['cases']), artifact_cases=2*len(result['cases']),
                   biological_truth_verified=False)
    (root/'results/passive_memory_verification.json').write_text(json.dumps(checked, indent=2)+'\n')
    print(json.dumps(checked, indent=2))


if __name__ == '__main__':
    verify()
