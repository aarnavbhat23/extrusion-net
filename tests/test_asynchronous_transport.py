import unittest
import numpy as np
from scipy.integrate import quad
from scipy.special import gammainc
from asynchronous_transport import reach_kernel,contacts,family_kernel,fit_local,fit_asymmetric
from study_load_probes import patterns


class AsynchronousTests(unittest.TestCase):
    def test_recurrence_matches_independent_continuous_time_quadrature(self):
        for mu in (.1,1.,3.):
            kernel=reach_kernel(5,mu)
            for a,b in ((0,0),(0,3),(1,1),(1,3),(3,2),(4,4)):
                def integrand(t):
                    pl=gammainc(a,t) if a else 1.;pr=gammainc(b,t) if b else 1.
                    return mu*np.exp(-mu*t)*pl*pr
                expected=quad(integrand,0,np.inf,epsabs=1e-12,epsrel=1e-12)[0]
                self.assertAlmostEqual(kernel[a,b],expected,places=11)

    def test_two_bands_are_indistinguishable_but_fourth_band_separates(self):
        rng=np.random.default_rng(121);m=10;s=.6;gates=rng.uniform(.2,.7,(2,m-1));loads=patterns(m,'alternating',None)
        maps=[contacts(q,*gates,family_kernel(m,s,'asynchronous')) for q in loads]
        second=np.stack([np.diag(c,2) for c in maps]);third=np.stack([np.diag(c,3) for c in maps])
        async_fit=fit_local(loads,second,third,'asynchronous');sync_fit=fit_local(loads,second,third,'synchronized')
        self.assertAlmostEqual(async_fit['survival'],s,places=12)
        self.assertAlmostEqual(sync_fit['survival'],(1+s)/2,places=12)
        for q,truth in zip(loads,maps):
            predicted=contacts(q,sync_fit['left'],sync_fit['right'],family_kernel(m,sync_fit['survival'],'synchronized'))
            np.testing.assert_allclose(np.diag(predicted,2),np.diag(truth,2),atol=1e-14)
            np.testing.assert_allclose(np.diag(predicted,3),np.diag(truth,3),atol=1e-14)
            self.assertGreater(np.max(np.diag(truth,4)-np.diag(predicted,4)),1e-5)

    def test_three_bands_through_distance_four_recover_unequal_rates_and_gates(self):
        rng=np.random.default_rng(122);m=14;loads=patterns(m,'alternating',None)
        for left_rate,right_rate,mu in ((1.,.5,.4),(1.,2.,1.),(.3,1.,.1)):
            gates=rng.uniform(.2,.7,(2,m-1));h=reach_kernel(m,mu,left_rate,right_rate)
            maps=[contacts(q,*gates,h) for q in loads]
            bands=[np.stack([np.diag(c,k) for c in maps]) for k in (2,3,4)]
            fit=fit_asymmetric(loads,*bands)
            self.assertAlmostEqual(fit['left_rate_over_removal'],left_rate/mu,places=9)
            self.assertAlmostEqual(fit['right_rate_over_removal'],right_rate/mu,places=9)
            np.testing.assert_allclose(np.stack([fit['left'],fit['right']]),gates,atol=1e-10)
            # Common scaling of all event rates cannot be observed in reach data.
            np.testing.assert_allclose(h,reach_kernel(m,mu*3,left_rate*3,right_rate*3),atol=1e-14)

    def test_fourth_band_formula_and_unequal_rates(self):
        for s in (.3,.8,.99):
            h=family_kernel(5,s,'asynchronous')
            self.assertAlmostEqual(h[1,1]/s**2,2/(1+s),places=12)
            self.assertAlmostEqual(h[1,2]/s**3,(s+3)/(1+s)**2,places=12)
        h=reach_kernel(5,.4,1.,2.)
        self.assertAlmostEqual(h[2,0],(1/1.4)**2)
        self.assertAlmostEqual(h[0,2],(2/2.4)**2)
        self.assertFalse(np.allclose(h,h.T))


if __name__=='__main__':unittest.main()
