import unittest

import torch

from models import extrusion_log_contacts
from polynomial_scan import (contact_argument_bound,contact_power_scan,
                             saturated_polynomial_scan,taylor_entry_error)


def dense_contacts(load,pl,pr):
    return extrusion_log_contacts(load.log(),pl.log(),pr.log()).exp()


class TestPolynomialScan(unittest.TestCase):
    def test_entrywise_powers_and_gradients_match_dense(self):
        torch.manual_seed(77)
        for degree in (1,2,4):
            args=[(torch.rand(2,9,dtype=torch.float64)*0.5+0.1).requires_grad_() for _ in range(3)]
            values=torch.randn(2,9,3,dtype=torch.float64,requires_grad=True)
            actual=contact_power_scan(*args,values,degree)
            expected=dense_contacts(*args)**degree@values
            torch.testing.assert_close(actual,expected,rtol=1e-10,atol=1e-12)
            probe=torch.randn_like(actual)
            actual_grads=torch.autograd.grad((actual*probe).sum(),args+[values],retain_graph=True)
            expected_grads=torch.autograd.grad((expected*probe).sum(),args+[values])
            for a,b in zip(actual_grads,expected_grads):
                torch.testing.assert_close(a,b,rtol=1e-9,atol=1e-11)

    def test_saturation_error_bound_and_hard_barriers(self):
        torch.manual_seed(8)
        load=torch.rand(13,dtype=torch.float64)*0.2
        pl=torch.rand(13,dtype=torch.float64)*0.7
        pr=torch.rand(13,dtype=torch.float64)*0.7
        values=torch.randn(13,3,dtype=torch.float64)
        beta=0.8
        dense=-torch.expm1(-beta*dense_contacts(load,pl,pr))
        bound=contact_argument_bound(load,pl,pr,beta)
        self.assertGreaterEqual(float(bound),float((beta*dense_contacts(load,pl,pr)).max()))
        for degree in (2,4,6):
            approximate=saturated_polynomial_scan(load,pl,pr,values,beta,degree)
            limit=taylor_entry_error(bound,degree)*values.abs().sum(0)
            self.assertTrue(torch.all((approximate-dense@values).abs()<=limit+1e-12))
        # The lifted scan, unlike a prefix quotient, permits exact zero passes.
        zeros=torch.zeros_like(pl)
        expected=torch.zeros(13,13,dtype=torch.float64)
        expected.diagonal(1).copy_(load[1:]**3)
        expected=expected+expected.T
        torch.testing.assert_close(contact_power_scan(load,zeros,zeros,values,3),expected@values)

    def test_large_argument_bound_is_not_hidden(self):
        ones=torch.ones(20,dtype=torch.float64)
        bound=contact_argument_bound(ones,ones,ones,4.)
        self.assertEqual(float(bound),76.)
        self.assertGreater(float(taylor_entry_error(bound,4)),1e6)


if __name__=='__main__':
    unittest.main()
