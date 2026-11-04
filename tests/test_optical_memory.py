import itertools
import unittest
import numpy as np
from study_optical_memory import contrast_interval, occupancy, observation_probability


class OpticalMemoryTests(unittest.TestCase):
    def test_independent_analytic_pulse_solution(self):
        for a, b in itertools.product([0., .2, 2.], [0., .4, 3.]):
            self.assertAlmostEqual(occupancy('AB', a, b), (1-np.exp(-a))*np.exp(-b))
            self.assertAlmostEqual(occupancy('BA', a, b), 1-np.exp(-a))

    def test_full_observation_law_by_explicit_hidden_paths(self):
        # Independently form the two analytic transition matrices, no expm.
        a, b = .6, 1.2
        matrices = {'A': np.array([[np.exp(-a), 1-np.exp(-a)], [0, 1]]),
                    'B': np.array([[1, 0], [1-np.exp(-b), np.exp(-b)]])}
        emission = [1., .4]
        for word in ['AB', 'BA', 'ABBA']:
            total = 0.
            for observed in itertools.product([0, 1], repeat=len(word)):
                probability = 0.
                for hidden in itertools.product([0, 1], repeat=len(word)):
                    weight, previous = 1., 0
                    for letter, state, value in zip(word, hidden, observed):
                        weight *= matrices[letter][previous, state]
                        weight *= emission[state] if value else 1-emission[state]
                        previous = state
                    probability += weight
                self.assertAlmostEqual(probability, observation_probability(word, observed, a, b, emission))
                total += probability
            self.assertAlmostEqual(total, 1.)

    def test_calibration_inverse_contains_every_admissible_contrast(self):
        for truth, sensitivity, offset in itertools.product(np.linspace(-1, 1, 21),
                                                            [.6, .7, .8], [-.02, 0., .02]):
            lower, upper = contrast_interval(sensitivity*truth+offset)
            self.assertLessEqual(lower, truth+1e-12)
            self.assertGreaterEqual(upper, truth-1e-12)
        with self.assertRaises(ValueError):
            contrast_interval(.2, calibration=(0., .8))

    def test_uncalibrated_crosstalk_is_exact_counterexample(self):
        # The biological and optical labels cannot alter an identical emission law.
        delta = occupancy('BA', .6, 1.2)-occupancy('AB', .6, 1.2)
        reporter_difference = .7*delta
        lower, _ = contrast_interval(reporter_difference)
        self.assertGreater(lower, 0)
        # Either biological occupancy OR an optical-state-sensitive reporter can
        # produce it; the latter violates the stated independent-reporter contract.
        self.assertAlmostEqual(reporter_difference, .7*(1-np.exp(-.6))*(1-np.exp(-1.2)))


if __name__ == '__main__':
    unittest.main()
