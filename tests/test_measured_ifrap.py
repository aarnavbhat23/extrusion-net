import unittest
import numpy as np
from study_measured_ifrap import fit, survival


class MeasuredIfrapTests(unittest.TestCase):
    def test_raw_signal_alias_and_sharp_intervention_endpoints(self):
        t = np.linspace(0, 122, 245)
        b1, k1, v1 = .004, .030, .008
        b2, k2, v2 = .009, .025, .003
        np.testing.assert_allclose(np.exp(-(b1+k1)*t), np.exp(-(b2+k2)*t), atol=1e-14)
        np.testing.assert_allclose(np.exp(-(b1+v1)*t), np.exp(-(b2+v2)*t), atol=1e-14)
        target = lambda b: np.exp(-(2*.034-.012-b)*20)
        candidates = [target(b) for b in np.linspace(0, .012, 1001)]
        self.assertAlmostEqual(min(candidates), np.exp(-1.12))
        self.assertAlmostEqual(max(candidates), np.exp(-.88))
        self.assertGreater(abs(target(b2)-target(b1)), .03)

    def test_mixture_reduces_to_single_and_plateau(self):
        t = np.arange(50, dtype=float)
        np.testing.assert_allclose(survival(t, [.03, .1, 1.], 'mixture'), np.exp(-.03*t))
        np.testing.assert_allclose(survival(t, [0., .03, .4], 'mixture'), .4+.6*np.exp(-.03*t))

    def test_observation_weighting_identity(self):
        rng = np.random.default_rng(515)
        x, z = rng.uniform(.1, 2, (2, 3, 40))
        prediction = np.linspace(1, .2, 40)
        np.testing.assert_allclose((x-z*prediction)**2, z*z*(x/z-prediction)**2, atol=1e-14)
        self.assertNotAlmostEqual(float(np.mean((x-z*prediction)**2)), float(np.mean((x/z-prediction)**2)))

    def test_known_noiseless_survival_with_variable_reference(self):
        t = np.linspace(0, 90, 80)
        expected = .35*np.exp(-.01*t)+.65*np.exp(-.08*t)
        z = np.stack([np.exp(-.003*t), np.exp(-.008*t)])
        for method in ['mixture', 'raw_contrast_mixture']:
            result = fit(t, z*expected, z, method)
            self.assertTrue(result['success'])
            self.assertEqual(len(result['starts']), 12)
            np.testing.assert_allclose(result['prediction'], expected, atol=2e-6)


if __name__ == '__main__':
    unittest.main()
