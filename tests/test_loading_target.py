import unittest
import numpy as np
from finite_loss import reference, TIMES
from loading_target import snapshot_bounds, timestamp_bounds, witness, binomial_interval, drift_bounds, drift_witness


class LoadingTargetTests(unittest.TestCase):
    def test_loading_drift_bounds_recover_endpoints_and_contain_true_target(self):
        ell, mu, delta, ratio = .8, .3, .1, .9
        means = [ell*.3/(.3+mu+delta), ratio*ell*.3/(.3+2*mu+delta)]
        result = drift_bounds([(p, p) for p in means])
        target = 1-ell*.3/(.3+3*mu+delta)
        self.assertLessEqual(result['interval'][0], target)
        self.assertGreaterEqual(result['interval'][1], target)
        for i, certificate in enumerate(result['certificates']):
            point = drift_witness(certificate['primal'])
            self.assertGreaterEqual(point['delta'], -1e-12)
            self.assertGreaterEqual(point['loading_ratio'], .8-1e-12)
            self.assertLessEqual(point['loading_ratio'], 1+1e-12)
            a = point['mu']+point['delta']
            self.assertAlmostEqual(.3/(.3+a), means[0])
            self.assertAlmostEqual(point['loading_ratio']*.3/(.3+a+(point['factor']-1)*point['mu']), means[1])
            self.assertAlmostEqual(point['target_absence'], result['interval'][i])

    def test_common_loading_is_a_special_case_of_drift_polytope(self):
        intervals = [(.3, .35), (.2, .25)]
        np.testing.assert_allclose(snapshot_bounds(intervals)['interval'], drift_bounds(intervals, ratio=(1., 1.))['interval'], atol=1e-12)

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

    def test_loading_gauge_preserves_full_snapshots_only_at_calibration_clock(self):
        ell, gamma = .7, .3
        p = np.array([1.1, .9, .25, .4, .7, 1., .1])
        scale = 1/ell
        other = p.copy()
        other[[0, 1, 2]] *= scale
        other[6] = scale*(p[6]+gamma)-gamma
        for factor in [1., 2., 3.]:
            first = ell*reference(tuple(p), factor, gamma)[:-1]
            second = reference(tuple(other), factor, gamma)[:-1]
            np.testing.assert_allclose(first, second, atol=1e-12)
        first = ell*reference(tuple(p), 3., .1)[:-1]
        second = reference(tuple(other), 3., .1)[:-1]
        expected = abs(ell*.1/(.1+3*p[2]+p[6])-.1/(.1+3*other[2]+other[6]))
        self.assertAlmostEqual(abs(first.sum()-second.sum()), expected, places=12)
        self.assertNotAlmostEqual(first.sum(), second.sum(), places=6)


if __name__ == '__main__':
    unittest.main()
