import math
import unittest
import numpy as np
from clock_calibration import flow_snapshot, generator_snapshot, odds_interval, speed_interval


class ClockCalibrationTests(unittest.TestCase):
    def test_independent_generator_matches_flow_at_boundaries(self):
        for p in ((1.2,.8,.3,.4,.6,.5), (.7,1.3,.8,0,1,1.2), (1,1,.4,0,0,.6)):
            for gamma in (.1, .3):
                np.testing.assert_allclose(flow_snapshot(p, gamma), generator_snapshot(p, gamma), atol=1e-13, rtol=0)

    def test_normalization_ambiguity_and_intervention_separation(self):
        p = (1.2,.8,.3,.4,.6,.5); gamma=.3; scale=1.7
        mu=scale*p[2]+(scale-1)*gamma
        changed=(scale*p[0],scale*p[1],mu,p[3],p[4],(scale*p[2]*p[5]+(scale-1)*gamma)/mu)
        a,b=flow_snapshot(p,gamma)[:-1],flow_snapshot(changed,gamma)[:-1]
        np.testing.assert_allclose(a/a.sum(),b/b.sum(),atol=1e-14,rtol=0)
        self.assertGreater(abs(flow_snapshot(p,gamma,2)[-1]-flow_snapshot(changed,gamma,2)[-1]),.01)

    def test_binomial_odds_interval_exact_coverage(self):
        for p in (.1,.4,.8):
            coverage=0.
            for x in range(16):
                lo,hi=odds_interval(x,15-x,.025)
                if lo <= p/(1-p) and (hi is None or p/(1-p) <= hi):
                    coverage += math.comb(15,x)*p**x*(1-p)**(15-x)
            self.assertGreaterEqual(coverage,.975-1e-14)

    def test_zero_counts_abstain_and_equal_clocks_fail(self):
        zero=np.zeros(81,int)
        self.assertEqual(speed_interval([(.1,zero),(.3,zero)])['bounds'],[0.,None])
        with self.assertRaises(ValueError):speed_interval([(.3,zero),(.3,zero)])

    def test_two_clock_exact_inverse_including_capture_boundaries(self):
        for parameters in ((1.15,.85,.27,.35,.65,.45),(.7,1.3,.8,0,1,1.2),
                           (1.3,.7,.4,1,0,.5),(1,1,.4,0,0,.6)):
            first = generator_snapshot(parameters,.1)[:-1].reshape(9,9)
            second = generator_snapshot(parameters,.3)[:-1].reshape(9,9)
            first = first/first.sum(); second = second/second.sum()
            a1,a2 = first[1,0]/first[0,0],second[1,0]/second[0,0]
            b1 = first[0,1]/first[0,0]
            vl = .2/(1/a2-1/a1); vr = vl*b1/a1
            mu = vl/a1-vl-vr-.1
            pl = 1-first[4,0]/(first[2,0]*a1*a1)
            pr = 1-first[0,4]/(first[0,2]*b1*b1)
            np.testing.assert_allclose((vl,vr,mu,pl,pr),parameters[:5],atol=1e-12,rtol=0)
            if parameters[3] > 0:
                d = first[3,0]/first[2,0]-(1-pl)*a1
                h = (vl*pl/d-vr-.1)/mu
            elif parameters[4] > 0:
                d = first[0,3]/first[0,2]-(1-pr)*b1
                h = (vr*pr/d-vl-.1)/mu
            else:
                # h has no effect when neither barrier captures.
                changed = (*parameters[:5],1.4)
                np.testing.assert_allclose(generator_snapshot(parameters,.3,2),generator_snapshot(changed,.3,2),atol=1e-13,rtol=0)
                continue
            self.assertAlmostEqual(h,parameters[5],places=11)

    def test_absolute_retention_recovers_speeds_with_one_clock(self):
        parameters = (1.15,.85,.27,.35,.65,.45); gamma = .3
        p = generator_snapshot(parameters,gamma)[:-1].reshape(9,9)
        total_rate = gamma/p[0,0]
        vl,vr = total_rate*p[1,0]/p[0,0],total_rate*p[0,1]/p[0,0]
        mu = total_rate-vl-vr-gamma
        np.testing.assert_allclose((vl,vr,mu),parameters[:3],atol=1e-12,rtol=0)


if __name__ == '__main__':unittest.main()
