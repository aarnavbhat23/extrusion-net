import unittest
import numpy as np
from finite_loss import HIGH, LOW, RATE, WEIGHTS, reference, scheduled, recover_unknown_factor
from study_loss_perturbations import snapshot


class FiniteLossTests(unittest.TestCase):
    def test_unknown_factor_recovery_from_independent_exact_snapshots(self):
        for h in (.4, .95, 1.05, 1.4):
            for pl, pr in ((.4, .7), (0., .7)):
                for factor in (.6, 1.8, 2.2):
                    p = (1.05, .95, .32, pl, pr, h, .12)
                    recovered = recover_unknown_factor(reference(p, 1., .3), reference(p, factor, .3))
                    np.testing.assert_allclose(recovered, (*p, factor), atol=1e-10, rtol=0)

    def test_equal_susceptibility_has_exact_intervention_ambiguity(self):
        first = (1.05, .95, .32, .4, .7, 1., .12)
        mu = .32*.8/1.2
        second = (1.05, .95, mu, .4, .7, 1., .44-mu)
        for gamma in (.1, .3, .8):
            np.testing.assert_allclose(reference(first, 1., gamma), reference(second, 1., gamma), atol=1e-12, rtol=0)
            np.testing.assert_allclose(reference(first, 1.8, gamma), reference(second, 2.2, gamma), atol=1e-12, rtol=0)
        span = abs(reference(first, 3., .3)[-1]-reference(second, 3., .3)[-1])
        self.assertGreater(span, .039)
        with self.assertRaises(ValueError):
            recover_unknown_factor(reference(first, 1., .3), reference(first, 1.8, .3))

    def test_uniformization_agrees_with_full_generator_at_domain_extremes(self):
        for p in (tuple(LOW), tuple(HIGH), (1.2, .8, .4, .6, .2, .7, .15)):
            for factor in (1., 2.):
                np.testing.assert_allclose(scheduled(p, factor), reference(p, factor), atol=2e-12, rtol=0)
        self.assertLess(abs(1 - WEIGHTS.sum(axis=0)).max(), 1e-13)
        self.assertLess(2 * HIGH[0] + 2 * HIGH[2] * max(1, HIGH[5]) + HIGH[6], RATE)

    def test_random_clock_flow_matches_full_resolvent(self):
        p = (1.2, .8, .4, .6, .2, .7, .15)
        for factor, gamma in ((1., .1), (1., .3), (2., .3), (3., .3)):
            np.testing.assert_allclose(snapshot(p, (1., factor, gamma), False), reference(p, factor, gamma), atol=1e-12, rtol=0)

    def test_state_dependent_loss_reduces_retention(self):
        p = (1.2, .8, .4, .6, .2, .7, .15)
        for gamma in (None, .3):
            baseline = np.atleast_2d(reference(p, 2., gamma))
            stressed = np.atleast_2d(reference(p, 2., gamma, .15))
            self.assertTrue(np.all(stressed[:, -1] > baseline[:, -1]))
            self.assertTrue(np.all(stressed[:, :-1] <= baseline[:, :-1] + 1e-12))

    def test_split_likelihood_coverage_by_exact_binomial_enumeration(self):
        # Any independent training estimator is allowed. Verify conditional
        # coverage for an intentionally bad fixed estimator, including zeros.
        from scipy.stats import binom
        for truth in (0., .1, .5, .9, 1.):
            for trained in (0., .01, .5, .99, 1.):
                probability = binom.pmf(np.arange(21), 20, truth)
                fitted = binom.pmf(np.arange(21), 20, trained)
                excluded_mass = probability[probability < .05 * fitted].sum()
                self.assertLessEqual(excluded_mass, .05 + 1e-12)


if __name__ == '__main__':
    unittest.main()
