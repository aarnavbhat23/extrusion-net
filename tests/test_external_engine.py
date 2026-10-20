import unittest
import numpy as np
from external_engine import binomial_interval,calibrate,predict


class ExternalBridgeTests(unittest.TestCase):
    def test_site_origin_bridge_and_noiseless_calibration(self):
        parameters=np.array([.7,.8,.3]);reach=predict(parameters);fit=calibrate(reach,10**12)
        self.assertFalse(fit['empty']);np.testing.assert_allclose(fit['point'],parameters,atol=1e-11)
        self.assertEqual(reach[0,0],1.)
        self.assertEqual(reach[4,0],.7**4*.8)
        self.assertTrue(np.all(predict(fit['low'])<=reach+1e-11))
        self.assertTrue(np.all(reach<=predict(fit['high'])+1e-11))
        self.assertTrue(np.all(predict(parameters,intervention='left')[4:]==0))

    def test_binomial_boundary_and_incompatible_calibration(self):
        self.assertEqual(binomial_interval(0,100,.05)[0],0.)
        self.assertEqual(binomial_interval(100,100,.05)[1],1.)
        zero=np.zeros((9,9));zero[0,0]=1.;fit=calibrate(zero,100)
        self.assertTrue(np.isfinite(predict(fit['low'])).all())
        self.assertTrue(np.isfinite(predict(fit['high'])).all())
        contradictory=np.ones((9,9));contradictory[1,0]=.2
        self.assertTrue(calibrate(contradictory,10000)['empty'])

    def test_unmodified_engine_episode_ends_before_reset(self):
        from external_engine import make_model,episodes
        reach,actions=episodes(make_model(.5,.2,.8),200,14000000)
        self.assertEqual(sum(actions.values()),200)
        self.assertEqual(actions['Evict'],0)
        self.assertEqual(reach[0,0],1.)
        self.assertTrue(np.all(np.diff(reach,axis=0)<=0))
        self.assertTrue(np.all(np.diff(reach,axis=1)<=0))


if __name__=='__main__':unittest.main()
