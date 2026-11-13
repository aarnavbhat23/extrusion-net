"""Check vector observations against the actual production scan and gradients."""
import unittest
import numpy as np
import torch
from scipy.optimize._numdiff import approx_derivative
from models import extrusion_linear_scan
from study_message_probes import message_response,shared_response
from study_load_probes import patterns


class MessageProbeTests(unittest.TestCase):
    def test_matches_production_scan(self):
        rng=np.random.default_rng(86);m=8
        q=rng.uniform(.1,.9,m);l=rng.uniform(.2,.8,m-1);r=rng.uniform(.2,.8,m-1)
        values=rng.normal(size=(2,m+1));output,jac=message_response(q,l,r,values)
        tensor=lambda a:torch.tensor(a,dtype=torch.float64)
        actual=extrusion_linear_scan(tensor(np.r_[.5,q]),tensor(np.r_[.5,l,.5]),
                                     tensor(np.r_[.5,r,.5]),tensor(values.T)).numpy().T.ravel()
        np.testing.assert_allclose(output,actual,atol=1e-13)
        numerical=approx_derivative(lambda g:message_response(q,g[:m-1],g[m-1:],values)[0],np.r_[l,r])
        np.testing.assert_allclose(jac,numerical,atol=1e-9,rtol=1e-6)

    def test_shared_gradient(self):
        rng=np.random.default_rng(87);x=rng.normal(size=(7,4));theta=rng.normal(size=8)
        q=patterns(8,'alternating',rng);values=rng.choice([-1.,1.],size=(2,9))
        _,jac=shared_response(theta,x,q,values)
        numerical=approx_derivative(lambda t:shared_response(t,x,q,values)[0],theta)
        np.testing.assert_allclose(jac,numerical,atol=1e-9,rtol=1e-6)


if __name__=='__main__':unittest.main()
