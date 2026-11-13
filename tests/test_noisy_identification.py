"""Independent numerical checks for noisy inverse fitting and local design."""
import unittest

import numpy as np

from study_identification import contact_array
from study_noisy_identification import contacts_and_jacobian, fit_map, local_design_scores, noise_aware_site


class NoisyIdentificationTests(unittest.TestCase):
    def test_forward_sensitivity_matches_independent_differences(self):
        rng=np.random.default_rng(77);size=6
        parameters=rng.uniform(.3,.8,3*size-2)
        values,jac=contacts_and_jacobian(parameters,size)
        expected=contact_array(parameters[:size],parameters[size:2*size-1],parameters[2*size-1:])
        np.testing.assert_allclose(values,expected[np.triu_indices(size)],atol=1e-14)
        for column in range(len(parameters)):
            plus=parameters.copy();minus=parameters.copy()
            plus[column]+=1e-6;minus[column]-=1e-6
            def independent(p):
                return contact_array(p[:size],p[size:2*size-1],p[2*size-1:])[np.triu_indices(size)]
            np.testing.assert_allclose(jac[:,column],(independent(plus)-independent(minus))/2e-6,atol=2e-9)

    def test_local_design_matches_direct_information_inverse(self):
        rng=np.random.default_rng(79);size=6
        parameters=rng.uniform(.4,.7,3*size-2)
        _,jac=contacts_and_jacobian(parameters,size)
        h=np.zeros((size-1,len(parameters)))
        for site in range(size-1):
            h[site,site+1]=parameters[size+site]
            h[site,size+site]=parameters[site+1]
        expected=[]
        for site in range(size-1):
            information=jac.T@jac+np.outer(h[site],h[site])
            variance=np.diag(h@np.linalg.solve(information,h.T))
            expected.append(np.mean(np.delete(variance,site)))
        np.testing.assert_allclose(local_design_scores(parameters,size),expected,rtol=1e-10)

    def test_exact_baseline_and_intervention_recover_heldout_responses(self):
        rng=np.random.default_rng(78);size=6
        parameters=rng.uniform(.4,.7,3*size-2)
        baseline=contact_array(parameters[:size],parameters[size:2*size-1],parameters[2*size-1:])
        value=parameters[1]*parameters[size]
        fit=fit_map(baseline,np.full_like(parameters,.55),(0,value))
        self.assertTrue(fit.success)
        np.testing.assert_allclose(fit.x,parameters,atol=2e-7)
        self.assertTrue(0<=noise_aware_site(parameters,size)<size-1)


if __name__=='__main__':unittest.main()
