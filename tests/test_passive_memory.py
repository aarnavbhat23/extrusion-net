import unittest
import numpy as np
from scipy.linalg import expm
from study_passive_memory import response, passive_envelope


class PassiveMemoryTests(unittest.TestCase):
    def test_equal_current_gaussian_law_can_have_different_passive_future(self):
        rates = np.array([.2, 2.])
        cov_a = np.eye(2); cov_b = np.eye(2); c = np.ones(2)
        means_a, means_b = np.array([1., -1.]), np.zeros(2)
        self.assertEqual(c @ means_a, c @ means_b)
        self.assertEqual(c @ cov_a @ c, c @ cov_b @ c)
        self.assertGreater(response(rates, means_a, [1.])[0], .5)
        self.assertEqual(response(rates, means_b, [1.])[0], 0.)

    def test_bound_holds_for_rotated_symmetric_dynamics(self):
        rotation = np.array([[1., -1.], [1., 1.]])/np.sqrt(2)
        a = rotation @ np.diag([.2, 3.]) @ rotation.T
        c = np.array([.4, 1.2]); initial_a = np.array([1., -.5]); initial_b = -initial_a
        times = np.linspace(0, 30, 101)
        observed = [abs(c @ expm(-a*t) @ (initial_a-initial_b)) for t in times]
        bound = passive_envelope(times, .2, np.linalg.norm(initial_a), np.linalg.norm(c))
        self.assertTrue(np.all(observed <= bound+1e-12))

    def test_two_mode_recurrence_is_defeated_by_another_passive_mode(self):
        times = np.arange(3.)
        a, b = np.exp(-.2), np.exp(-2.)
        contrast = [a*b, -(a+b), 1.]
        null = response([.2, 2.], [1., -1.], times)
        self.assertAlmostEqual(np.dot(contrast, null), 0.)
        artifact = response([.2, 2., .7], [1., -1., .4], times)
        self.assertGreater(abs(np.dot(contrast, artifact)), .01)


if __name__ == '__main__':
    unittest.main()
