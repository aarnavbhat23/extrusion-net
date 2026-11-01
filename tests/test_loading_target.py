import unittest
import numpy as np
from finite_loss import reference, TIMES
from loading_target import snapshot_bounds, timestamp_bounds, witness, binomial_interval


class LoadingTargetTests(unittest.TestCase):
    def test_loading_cancellation_agrees_with_full_generator(self):
        for loading in [.6, .8, 1.]:
            p = (1.1, .9, .3, .4, .7, 1., .15)
            r0 = loading*reference(p, 1., .3)[:-1].sum()
            r1 = loading*reference(p, 2., .3)[:-1].sum()
            target = loading*reference(p, 3., .3)[:-1].sum()
            self.assertAlmostEqual(1/target, 1/r0+2*(1/r1-1/r0), places=11)
            for f in [1., 2., 3.]:
                observed = loading*reference(p, f)[:, :-1].sum(axis=1)
                np.testing.assert_allclose(observed, loading*np.exp(-(f*.3+.15)*np.array(TIMES)), atol=1e-12)

    def test_exact_data_sharp_endpoints_have_feasible_mechanisms(self):
        r0, r1 = .8*.3/(.3+.4), .8*.3/(.3+.7)
        result = snapshot_bounds([(r0, r0), (r1, r1)])
        endpoints = []
        for certificate in result['certificates']:
            mechanism = witness(certificate['primal'])
            self.assertGreaterEqual(mechanism['delta'], -1e-12)
            self.assertGreaterEqual(mechanism['factor'], 1.8-1e-12)
            self.assertLessEqual(mechanism['factor'], 2.2+1e-12)
            a = mechanism['mu']+mechanism['delta']
            self.assertAlmostEqual(.3/(.3+a), r0)
            self.assertAlmostEqual(.3/(.3+a+(mechanism['factor']-1)*mechanism['mu']), r1)
            endpoints.append(mechanism['target_absence'])
        np.testing.assert_allclose(endpoints, result['interval'], atol=1e-12)

    def test_sharp_snapshot_polytope_contains_dense_direct_mechanisms(self):
        intervals = [(.25, .4), (.15, .3)]
        result = snapshot_bounds(intervals)
        count = 0
        for ell in [.6, .8, 1.]:
            for mu in np.linspace(0.05, .8, 9):
                for delta in np.linspace(0, .4, 7):
                    for f in [1.8, 2., 2.2]:
                        r0 = ell*.3/(.3+mu+delta); r1 = ell*.3/(.3+f*mu+delta)
                        if not (intervals[0][0] <= r0 <= intervals[0][1] and intervals[1][0] <= r1 <= intervals[1][1]):
                            continue
                        q = 1-ell*.3/(.3+3*mu+delta)
                        self.assertLessEqual(result['interval'][0]-1e-12, q)
                        self.assertGreaterEqual(result['interval'][1]+1e-12, q)
                        count += 1
        self.assertGreater(count, 100)

    def test_timestamp_outer_interval_encloses_exact_target(self):
        ell, a, mu, f = .8, .4, .3, 2.
        means = [ell*np.exp(-rate*t) for rate in [a, a+(f-1)*mu] for t in [1/.3, 2/.3]]
        result = timestamp_bounds([(p, p) for p in means])
        target = 1-ell*.3/(.3+a+2*mu)
        self.assertLessEqual(result['interval'][0], target+1e-12)
        self.assertGreaterEqual(result['interval'][1], target-1e-12)

    def test_zero_counts_are_not_replaced_by_finite_reciprocal_cutoffs(self):
        interval = binomial_interval(0, 100, .025)
        result = snapshot_bounds([interval, (.1, .2)])
        self.assertEqual(result['interval'], [0., 1.])


if __name__ == '__main__':
    unittest.main()
