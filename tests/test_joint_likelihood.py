import itertools
import math
import unittest
import numpy as np
from joint_likelihood import multinomial_log_likelihood, split_counts


class JointLikelihoodTests(unittest.TestCase):
    def test_zero_events_are_not_clipped(self):
        self.assertEqual(multinomial_log_likelihood([0, 3], [0, 1]), 0)
        self.assertEqual(multinomial_log_likelihood([1, 2], [0, 1]), -np.inf)
        with self.assertRaises(ValueError):
            multinomial_log_likelihood([1, 2], [.2, .7])

    def test_split_conserves_counts_and_replays(self):
        counts = np.array([0, 5, 20, 9])
        a, b = split_counts(counts, np.random.default_rng(33))
        np.testing.assert_array_equal(a + b, counts)
        np.testing.assert_array_equal(a, split_counts(counts, np.random.default_rng(33))[0])
        self.assertEqual(a.sum(), counts.sum() // 2)

    def test_exact_split_likelihood_coverage_on_all_binomial_outcomes(self):
        # Enumerate training/test counts, including boundary fits and zero support.
        alpha = .05
        for p in (.05, .3, .5, .9):
            rejected = 0.
            for x, y in itertools.product(range(4), range(6)):
                q = x / 3
                probability = math.comb(3, x) * p**x * (1-p)**(3-x)
                probability *= math.comb(5, y) * p**y * (1-p)**(5-y)
                lp = multinomial_log_likelihood([y, 5-y], [p, 1-p])
                lq = multinomial_log_likelihood([y, 5-y], [q, 1-q])
                rejected += probability * (lp < lq + math.log(alpha))
            self.assertLessEqual(rejected, alpha + 1e-14)

    def test_ratio_expectation_with_zero_support_is_at_most_one(self):
        for p, q in ((.3, .8), (1., .5), (0., .8)):
            expectation = 0.
            for x in range(4):
                probability = math.comb(3, x) * p**x * (1-p)**(3-x)
                if probability == 0:
                    continue
                lp = multinomial_log_likelihood([x, 3-x], [p, 1-p])
                lq = multinomial_log_likelihood([x, 3-x], [q, 1-q])
                expectation += probability * math.exp(lq-lp)
            self.assertLessEqual(expectation, 1 + 1e-14)


if __name__ == '__main__':
    unittest.main()
