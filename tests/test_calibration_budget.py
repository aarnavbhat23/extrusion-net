"""Check the balanced-routing information bound against independent contacts."""
import unittest

import numpy as np

from study_calibration_budget import diffusion_bound, information_covariance, laplacian, tangent_contacts
from study_identification import contact_array


class CalibrationBudgetTests(unittest.TestCase):
    def test_tangent_matches_nonlinear_contact_derivatives(self):
        gates=7;a=.6;jac=tangent_contacts(gates,a)
        positions=[(i,i+d) for d in range(2,gates+1) for i in range(gates+1-d)]
        for k in range(gates):
            z=np.zeros(gates);z[k]=1e-6
            plus=contact_array(np.ones(gates+1),a+z,a-z)
            minus=contact_array(np.ones(gates+1),a-z,a+z)
            derivative=np.array([(plus[i,j]-minus[i,j])/2e-6 for i,j in positions])
            np.testing.assert_allclose(jac[:,k],derivative,atol=2e-9)

    def test_series_and_spectral_order(self):
        for a in (.3,.6,.8):
            series=sum(a**(2*d-2)*(d*(d*d-1)/6)**2 for d in range(2,1000))
            self.assertAlmostEqual(diffusion_bound(a)/series,1.,places=12)
            jac=tangent_contacts(16,a);fisher=jac.T@jac;path=laplacian(16)
            self.assertGreaterEqual(np.linalg.eigvalsh(fisher-a*a*path).min(),-1e-10)
            self.assertGreaterEqual(np.linalg.eigvalsh(diffusion_bound(a)*path-fisher).min(),-1e-8)
            np.testing.assert_allclose(jac@np.ones(16),0.,atol=1e-14)

    def test_anchored_path_inverse_and_variance_bounds(self):
        gates=12;a=.6;c=diffusion_bound(a)
        info=c*laplacian(gates);info[0,0]+=1
        np.testing.assert_allclose(np.diag(np.linalg.inv(info)),1+np.arange(gates)/c,atol=1e-12)
        covariance=information_covariance(tangent_contacts(gates,a),[0])
        variance=np.trace(covariance)/gates
        self.assertGreaterEqual(variance,1+(gates-1)/(2*c)-1e-12)
        self.assertLessEqual(variance,1+(gates-1)/(2*a*a)+1e-12)


if __name__=='__main__':unittest.main()
