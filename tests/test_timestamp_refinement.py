import unittest
import numpy as np
from study_timestamp_refinement import refine


class TimestampRefinementTests(unittest.TestCase):
    def test_slicing_preserves_analytic_extrema_and_controls_relaxation(self):
        # A=L+1, D=W=0, L in [-.5,0]. Presence exp(L)/(L+2) is increasing.
        original = dict(status='ok', problem=dict(matrix=[[-1., 1., 0., 0.], [1., -1., 0., 0.]],
                        rhs=[1., -1.], bounds=[[-.5, 0.], [0., None], [0., 0.], [0., 0.]]),
                        certificates=[dict(value=-.5), dict(value=0.)])
        result = refine(original)
        truth = [1-.5, 1-np.exp(-.5)/1.5]
        self.assertLessEqual(result['interval'][0], truth[0]+1e-12)
        self.assertGreaterEqual(result['interval'][1], truth[1]-1e-12)
        self.assertLessEqual(max(abs(np.array(result['interval'])-truth)), result['one_sided_relaxation_slack_bound']+1e-12)
        self.assertEqual(len(result['bins']), 16)


if __name__ == '__main__':
    unittest.main()
