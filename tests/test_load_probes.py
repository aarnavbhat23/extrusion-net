"""Independent inversion, sensitivity and analytic-gradient checks."""
import unittest
import numpy as np
from scipy.optimize._numdiff import approx_derivative
from study_identification import contact_array
from study_load_probes import patterns,recover_local,controller_output


class LoadProbeTests(unittest.TestCase):
    def test_two_views_identify_original_contact_sum(self):
        rng=np.random.default_rng(25)
        for size in (3,8,32):
            q=patterns(size,'alternating',rng)
            l=rng.uniform(.1,.9,size-1);r=rng.uniform(.1,.9,size-1)
            observed=np.stack([np.diag(contact_array(row,l,r),1) for row in q])
            a,b=recover_local(q,observed)
            np.testing.assert_allclose(a,l,atol=1e-14)
            np.testing.assert_allclose(b,r,atol=1e-14)

    def test_controller_gradient(self):
        rng=np.random.default_rng(26);x=rng.normal(size=(6,4));theta=rng.normal(size=8)
        q=patterns(7,'alternating',rng)
        for full in (False,True):
            _,jac=controller_output(theta,x,q,full)
            numerical=approx_derivative(lambda t:controller_output(t,x,q,full)[0],theta)
            np.testing.assert_allclose(jac,numerical,rtol=1e-6,atol=1e-9)

    def test_view_drift_bound_and_uniform_failure(self):
        rng=np.random.default_rng(27);q=patterns(16,'alternating',rng);epsilon=.02
        l=rng.uniform(.2,.8,15);r=rng.uniform(.2,.8,15)
        y=q[:,1:]*(l+rng.uniform(-epsilon,epsilon,(2,15)))+q[:,:-1]*(r+rng.uniform(-epsilon,epsilon,(2,15)))
        a,b=recover_local(q,y)
        self.assertLessEqual(np.max(abs(np.r_[a-l,b-r])),epsilon/.8+1e-14)
        uniform=patterns(16,'uniform',rng)
        y=uniform[:,1:]*l+uniform[:,:-1]*r
        a,b=recover_local(uniform,y)
        np.testing.assert_allclose(a,b,atol=1e-14)
        self.assertGreater(np.max(abs(a-l)),.01)


if __name__=='__main__':unittest.main()
