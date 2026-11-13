"""Independent production-scan and numerical-extremum checks."""
import unittest

import numpy as np
import torch
from scipy.optimize import minimize_scalar

from counterfactual_bounds import (cut_endpoints, initial_for_mix, mixer_outputs,
                                  sharp_output_bounds, symmetric_contacts)
from models import extrusion_linear_scan
from study_identification import contact_array, factor_at_initial, observable_bands


class CounterfactualOutputTests(unittest.TestCase):
    def test_segment_and_output_match_production_scan(self):
        rng=np.random.default_rng(89);size=7;site=3
        q=rng.uniform(.4,.9,size);left,right=rng.uniform(.4,.7,(2,size-1))
        baseline=contact_array(q,left,right)
        endpoints,interval,prefix=cut_endpoints(baseline,site,.3)
        noisy=baseline.copy();noisy[0,-1]+=.01
        with self.assertRaises(ValueError):cut_endpoints(noisy,site,.3)
        q,sums,products=observable_bands(baseline)
        values=rng.normal(size=(size+1,4));weight=rng.normal(size=(3,5));bias=rng.normal(size=3)
        bounds=sharp_output_bounds(endpoints,values,weight,bias)
        for theta in np.linspace(0,1,17):
            initial=initial_for_mix(theta,interval,prefix)
            l,r=factor_at_initial(q,sums,products,initial);r[site]*=.3
            actual=symmetric_contacts(contact_array(q,l,r))
            np.testing.assert_allclose(actual,(1-theta)*endpoints[0]+theta*endpoints[1],atol=1e-12)
            # Explicit mapping to the physical positions used by the production scan.
            load=torch.tensor(np.r_[.5,q],dtype=torch.float64)
            pl=torch.tensor(np.r_[.5,l,.5],dtype=torch.float64)
            pr=torch.tensor(np.r_[.5,r,.5],dtype=torch.float64)
            v=torch.tensor(np.c_[values,np.ones(size+1)],dtype=torch.float64)
            mixed=extrusion_linear_scan(load,pl,pr,v)
            features=torch.cat([mixed[:,:-1],torch.log1p(mixed[:,-1:])],dim=1)
            output=(features@torch.tensor(weight.T)+torch.tensor(bias)).numpy()
            np.testing.assert_allclose(output,mixer_outputs(actual,values,weight,bias),atol=1e-12)
            self.assertTrue(np.all(output>=bounds['lower']-1e-12))
            self.assertTrue(np.all(output<=bounds['upper']+1e-12))

    def test_interior_extremum_and_independent_scalar_optimization(self):
        endpoints=np.array([[[0.,0.],[0.,0.]],[[0.,1.],[1.,0.]]])
        values=np.array([[-1.],[-1.]])
        weight=np.array([[1.,1.5],[-1.,-1.5]])
        bounds=sharp_output_bounds(endpoints,values,weight)
        self.assertAlmostEqual(bounds['upper_theta'][0,0],.5)
        self.assertGreater(bounds['upper'][0,0],max(mixer_outputs(c,values,weight)[0,0] for c in endpoints)+.05)
        for channel in range(2):
            def output(theta):return mixer_outputs((1-theta)*endpoints[0]+theta*endpoints[1],values,weight)[0,channel]
            low=min(output(0),output(1),minimize_scalar(output,bounds=(0,1),method='bounded').fun)
            high=max(output(0),output(1),-minimize_scalar(lambda x:-output(x),bounds=(0,1),method='bounded').fun)
            self.assertAlmostEqual(bounds['lower'][0,channel],low,places=9)
            self.assertAlmostEqual(bounds['upper'][0,channel],high,places=9)
        same=sharp_output_bounds(np.stack([endpoints[1],endpoints[1]]),values,weight)
        np.testing.assert_allclose(same['lower'],same['upper'])


if __name__=='__main__':unittest.main()
