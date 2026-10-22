"""Independent checks of observational equivalence and intervention calibration."""
import unittest

import numpy as np
import torch

from screen_directions import contacts,equivalent
from study_identification import admissible_initial_interval,contact_array,factor_at_initial,gate_family,observable_bands


class IdentificationTests(unittest.TestCase):
    def test_distinct_gates_same_map_different_intervention(self):
        q=torch.tensor([.6,.7,.8,.65],dtype=torch.float64)
        left=torch.tensor([.55,.6,.5],dtype=torch.float64);right=torch.tensor([.6,.5,.65],dtype=torch.float64)
        new_l,new_r=equivalent(q,left,right,.7)
        self.assertTrue(bool(((new_l>0)&(new_l<1)&(new_r>0)&(new_r<1)).all()))
        baseline=contacts(q,left,right)
        torch.testing.assert_close(baseline,contacts(q,new_l,new_r),atol=1e-14,rtol=1e-14)
        rcut=right.clone();rcut[0]=0
        new_cut=new_r.clone();new_cut[0]=0
        self.assertGreater(float((contacts(q,left,rcut)-contacts(q,new_l,new_cut)).abs().max()),.01)
        np.testing.assert_allclose(baseline.numpy(),contact_array(q.numpy(),left.numpy(),right.numpy()),atol=1e-14)

    def test_inverse_bands_and_single_measurement(self):
        q=torch.tensor([.6,.7,.8,.65],dtype=torch.float64)
        left=torch.tensor([.55,.6,.5],dtype=torch.float64);right=torch.tensor([.6,.5,.65],dtype=torch.float64)
        c=contacts(q,left,right);inv=torch.linalg.inv(c)
        torch.testing.assert_close(torch.diag(inv,1),-left/q[:-1]-right/q[1:])
        torch.testing.assert_close(torch.diag(inv,2),right[:-1]*left[1:]/q[1:-1])
        self.assertLess(float(inv[0,3].abs()),1e-14)
        rcut=right.clone();rcut[0]=0
        recovered=(c[0,1]-contacts(q,left,rcut)[0,1])/q[0]
        l,r=equivalent(q,left,right,recovered)
        torch.testing.assert_close(l,left);torch.testing.assert_close(r,right)

    def test_feasible_interval_and_sharp_directional_envelopes(self):
        rng=np.random.default_rng(45)
        for size in (4,8,16):
            q=rng.uniform(.4,.9,size);left=rng.uniform(.4,.7,size-1);right=rng.uniform(.4,.7,size-1)
            c=contact_array(q,left,right)
            observed_q,sums,products=observable_bands(c)
            low,high=admissible_initial_interval(observed_q,sums,products)
            self.assertLessEqual(low,right[0]);self.assertGreaterEqual(high,right[0])
            endpoints=[factor_at_initial(q,sums,products,x) for x in (low,high)]
            for l,r in endpoints:
                self.assertTrue(np.min(l)>-1e-12 and np.max(l)<1+1e-12)
                self.assertTrue(np.min(r)>-1e-12 and np.max(r)<1+1e-12)
                np.testing.assert_allclose(contact_array(q,l,r),c,atol=1e-12)
            for site in range(size-1):
                boundary=[]
                for l,r in endpoints:
                    cut=r.copy();cut[site]=0
                    boundary.append(contact_array(q,l,cut))
                for initial in np.linspace(low,high,21):
                    l,r=factor_at_initial(q,sums,products,initial)
                    r[site]=0
                    intervention=contact_array(q,l,r)
                    self.assertTrue(np.all(intervention<=boundary[0]+1e-12))
                    self.assertTrue(np.all(intervention>=boundary[1]-1e-12))
            # Both-arm gate scalings have the same response in every factorization.
            multiplier=rng.uniform(.1,.9,size-1)
            changes=[contact_array(q,l*multiplier,r*multiplier) for l,r in endpoints]
            np.testing.assert_allclose(changes[0],changes[1],atol=1e-12)


if __name__=='__main__':unittest.main()
