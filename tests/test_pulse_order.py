import unittest
import numpy as np
from scipy.stats import fisher_exact
from pulse_order import (schedule,moments,mixture_survival,generators,endpoint_survival,
                         equal_arm_fisher,robust_threshold,robust_fisher_p)

class PulseOrderTests(unittest.TestCase):
    def test_exact_polynomial_moment_cancellation(self):
        for degree in range(5):
            w=schedule(degree)
            np.testing.assert_allclose(moments(w,degree),moments(1-w,degree),rtol=1e-13,atol=1e-12)

    def test_arbitrary_correlated_frailty_cancels(self):
        rng=np.random.default_rng(4)
        for degree in range(4):
            c=rng.uniform(0,.1,(13,2,degree+1));p=rng.dirichlet(np.ones(13));w=schedule(degree)
            self.assertAlmostEqual(mixture_survival(w,c,p),mixture_survival(1-w,c,p),places=13)

    def test_degree_beyond_null_can_mimic_memory(self):
        c=np.array([[[.02,0,.05],[.02,0,0]]]);w=schedule(1)
        self.assertGreater(abs(mixture_survival(w,c,[1])-mixture_survival(1-w,c,[1])),.01)
        w=schedule(2);self.assertAlmostEqual(mixture_survival(w,c,[1]),mixture_survival(1-w,c,[1]),places=13)

    def test_hidden_state_may_be_unobservable(self):
        q=generators(.5,.1,equal_hazards=True)
        self.assertGreater(np.linalg.norm(q[0]@q[1]-q[1]@q[0]),.01)
        for degree in range(3):
            w=schedule(degree);self.assertAlmostEqual(endpoint_survival(w,q),endpoint_survival(1-w,q),places=13)

    def test_equal_arm_fisher_matches_reference(self):
        for n in [5,30]:
            for x in range(n+1):
                for y in range(n+1):
                    expected=fisher_exact([[x,n-x],[y,n-y]]).pvalue
                    self.assertAlmostEqual(float(equal_arm_fisher(x,y,n)),expected,places=12)

    def test_supremum_fisher_matches_every_compatible_count(self):
        n=6
        for errors in (0,1,2):
            for x in range(n+1):
                for y in range(n+1):
                    values=[fisher_exact([[a,n-a],[b,n-b]]).pvalue
                        for a in range(max(0,x-errors),min(n,x+errors)+1)
                        for b in range(max(0,y-errors),min(n,y+errors)+1)]
                    self.assertAlmostEqual(float(robust_fisher_p(x,y,n,errors)),max(values),places=12)

    def test_corruption_inflation_and_empty_input(self):
        self.assertAlmostEqual(robust_threshold(1000,1000,5,5)-robust_threshold(1000,1000),.01)
        with self.assertRaises(ValueError):robust_threshold(0,10)
        with self.assertRaises(ValueError):schedule(-1)

if __name__=='__main__':unittest.main()
