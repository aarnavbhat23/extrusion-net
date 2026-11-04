"""Independent scan checks, alias sums and global fixed-point bounds."""
import unittest
import itertools
import numpy as np
import torch
from models import extrusion_linear_scan
from colored_recovery import (alias_contraction,choose_colors,colored_response,
                              local_inverse,recover,scan_messages)
from study_probe_mismatch import positive_message_bounds
from study_identification import contact_array
from counterfactual_bounds import symmetric_contacts


class ColoredRecoveryTests(unittest.TestCase):
    def test_scan_matches_production_for_arbitrary_signed_values(self):
        rng=np.random.default_rng(119);m=13;q=rng.uniform(.1,.9,m)
        l=rng.uniform(0,.8,m-1);r=rng.uniform(0,.8,m-1);v=rng.normal(size=(m+1,7))
        tensor=lambda x:torch.tensor(x,dtype=torch.float64)
        actual=extrusion_linear_scan(tensor(np.r_[0,q]),tensor(np.r_[0,l,0]),tensor(np.r_[0,r,0]),tensor(v)).numpy()
        np.testing.assert_allclose(scan_messages(q,l,r,v),actual,atol=1e-13)

    def test_alias_series_matches_direct_sum(self):
        for rho,k in ((.2,8),(.5,16),(.8,44),(.95,240)):
            direct=sum((t*k-2)*(t*k-3)*rho**(t*k-4)+(t*k+2)*(t*k+1)*rho**(t*k) for t in range(1,100))
            self.assertAlmostEqual(alias_contraction(rho,k),1.8/1.6*direct,places=12)

    def test_contraction_and_recovery_including_noise(self):
        rng=np.random.default_rng(120)
        for size,rho in ((12,.8),(64,.5),(128,.8)):
            k,rate=choose_colors(rho,size);truth=rng.uniform(0,rho,2*(size-1));other=rng.uniform(0,rho,2*(size-1))
            a=local_inverse(colored_response(truth,size,k),size,k)-truth
            b=local_inverse(colored_response(other,size,k),size,k)-other
            self.assertLessEqual(np.max(abs(a-b)),rate*np.max(abs(truth-other))+1e-12)
            clean=colored_response(truth,size,k);noise=rng.uniform(-.001,.001,clean.shape)
            estimate,fit=recover(clean+noise,size,rho,k,rate)
            bound=fit['fixed_point_error_bound']+.001/(.8*(1-rate))
            self.assertTrue(fit['converged']);self.assertLessEqual(np.max(abs(estimate-truth)),bound+1e-12)
            exact,result=recover(clean,size,rho,k,rate)
            self.assertLess(np.max(abs(exact-truth)),1e-9)
            self.assertGreaterEqual(result['fixed_point_error_bound']+1e-12,np.max(abs(exact-truth)))

    def test_positive_response_bounds_attained_at_gate_box_corners(self):
        rng=np.random.default_rng(121);size=4;center=rng.uniform(.1,.4,6);q=rng.uniform(.1,.9,4)
        low,high=positive_message_bounds(q,center,size,.05,.5);outputs=[]
        for bits in itertools.product((-1,1),repeat=6):
            g=center+.05*np.array(bits)
            outputs.append(symmetric_contacts(contact_array(q,g[:3],g[3:])).sum(axis=1))
        np.testing.assert_allclose(low,np.min(outputs,axis=0),atol=1e-13)
        np.testing.assert_allclose(high,np.max(outputs,axis=0),atol=1e-13)


if __name__=='__main__':unittest.main()
