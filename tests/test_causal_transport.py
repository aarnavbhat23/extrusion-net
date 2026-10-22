import unittest
import numpy as np
import torch
from causal_transport import causal_sum,clock_messages,CausalModel
from study_kinetic_transfer import kinetic_contacts


class CausalTransportTests(unittest.TestCase):
    def test_window_scan_and_kinetic_messages_match_independent_generator(self):
        torch.manual_seed(71);m=10
        load=torch.rand(m+1,dtype=torch.float64);left=torch.rand(m+1,dtype=torch.float64);right=torch.rand(m+1,dtype=torch.float64)
        v=torch.randn(m+1,3,dtype=torch.float64)
        for shared in (False,True):
            expected=np.tril(kinetic_contacts(load.numpy()[1:],left.numpy()[1:-1],right.numpy()[1:-1],.7,shared))@v.numpy()
            got=clock_messages(load,left,right,v,torch.tensor(.7,dtype=torch.float64),m,shared)
            np.testing.assert_allclose(got.detach().numpy(),expected,atol=2e-14)
        for width in (1,2,3,6,20):
            got=causal_sum(v,left,width);expected=torch.zeros_like(v)
            for i in range(m+1):
                for j in range(max(0,i-width+1),i+1):expected[i]+=v[j]*torch.prod(left[j:i])
            torch.testing.assert_close(got,expected)

    def test_prefix_invariance_including_local_backbone(self):
        for variant in ('static','shared','independent','attention'):
            model=CausalModel(variant,25,8).eval();x=torch.randint(25,(2,25));changed=x.clone();changed[:,13:]=torch.randint(25,(2,12))
            with torch.no_grad():torch.testing.assert_close(model(x)[:,:13],model(changed)[:,:13])

    def test_clock_gradients_match_finite_differences(self):
        torch.manual_seed(72)
        values=torch.randn(6,2,dtype=torch.float64)
        args=[torch.rand(6,dtype=torch.float64,requires_grad=True) for _ in range(3)]
        s=torch.tensor(.65,dtype=torch.float64,requires_grad=True)
        self.assertTrue(torch.autograd.gradcheck(lambda a,l,r,p:clock_messages(a,l,r,values,p,3,True),(*args,s)))


if __name__=='__main__':unittest.main()
