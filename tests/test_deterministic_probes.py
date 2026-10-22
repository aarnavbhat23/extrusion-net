"""Zero-gate isometry, nonlocal secant bounds and production observations."""
import unittest
import numpy as np
import torch
from models import extrusion_linear_scan
from deterministic_probes import coded_values,probe_response,recovery_constants,recovery_error_bound
from study_load_probes import patterns


class DeterministicProbeTests(unittest.TestCase):
    def test_zero_gate_singular_values_are_length_independent(self):
        for size in (2,3,8,17):
            _,j=probe_response(np.zeros(2*(size-1)),size)
            singular=np.linalg.svd(j,compute_uv=False)
            np.testing.assert_allclose(singular[:size-1],2.,atol=1e-13)
            np.testing.assert_allclose(singular[size-1:],1.6,atol=1e-13)

    def test_uniform_tail_and_global_secant_bounds(self):
        rng=np.random.default_rng(110)
        for size in (4,9,24):
            _,j0=probe_response(np.zeros(2*(size-1)),size)
            for rho in (.02,.05,.08):
                c=recovery_constants(rho);self.assertGreater(c['lower'],0)
                a=rng.uniform(0,rho,2*(size-1));b=rng.uniform(0,rho,2*(size-1))
                y,j=probe_response(a,size);z,_=probe_response(b,size)
                self.assertLessEqual(np.linalg.norm(j-j0,2),c['tail']+1e-12)
                ratio=np.linalg.norm(y-z)/np.linalg.norm(a-b)
                self.assertGreaterEqual(ratio,c['lower']-1e-12)
                self.assertLessEqual(ratio,c['upper']+1e-12)
                noise=rng.normal(0,.01,len(y));observation=z+noise
                bound=recovery_error_bound(np.linalg.norm(y-observation),np.linalg.norm(noise),rho)
                self.assertLessEqual(np.linalg.norm(a-b),bound+1e-12)
        with self.assertRaises(ValueError):recovery_error_bound(.1,.1,.5)

    def test_coded_responses_match_production_scan(self):
        rng=np.random.default_rng(111);size=10;gates=rng.uniform(0,.05,2*(size-1))
        values=coded_values(size);actual=[];tensor=lambda x:torch.tensor(x,dtype=torch.float64)
        for q in patterns(size,'alternating',None):
            result=extrusion_linear_scan(tensor(np.r_[.5,q]),tensor(np.r_[.5,gates[:size-1],.5]),
                                         tensor(np.r_[.5,gates[size-1:],.5]),tensor(values.T))
            actual.append(result.numpy().T.ravel())
        np.testing.assert_allclose(probe_response(gates,size)[0],np.concatenate(actual),atol=1e-13)


if __name__=='__main__':unittest.main()
