"""Enumeration controls and a coupled-clock band identity."""
import unittest
import numpy as np
from study_kinetic_transfer import kinetic_contacts,shared_clock_signature
from study_load_probes import patterns
from study_identification import contact_array
from counterfactual_bounds import symmetric_contacts


class KineticTests(unittest.TestCase):
    def test_independent_lifetimes_equal_scaled_gate_operator(self):
        rng=np.random.default_rng(130);q=rng.uniform(.1,.9,8);l=rng.uniform(.1,.7,7);r=rng.uniform(.1,.7,7)
        for s in (.5,.8,1.):
            expected=symmetric_contacts(contact_array(q,s*l,s*r))
            np.testing.assert_allclose(kinetic_contacts(q,l,r,s,False),expected,atol=1e-13)
        np.testing.assert_allclose(kinetic_contacts(q,l,r,1.,True),expected,atol=1e-13)

    def test_shared_lifetime_has_maximum_time_factor_and_signature(self):
        rng=np.random.default_rng(131);m=9;l=rng.uniform(.2,.7,m-1);r=rng.uniform(.2,.7,m-1)
        loads=patterns(m,'alternating',None)
        for s in (.5,.8,.95):
            maps=[kinetic_contacts(q,l,r,s,True) for q in loads]
            # Independent direct reach formula, not lifetime-mixture enumeration.
            expected=np.zeros((m+1,m+1));q=loads[0]
            for i in range(m):
                for j in range(i+1,m+1):
                    for bond in range(i,j):
                        expected[i,j]+=q[bond]*np.prod(l[i:bond])*np.prod(r[bond:j-1])*s**max(bond-i,j-1-bond)
            np.testing.assert_allclose(maps[0],expected+expected.T,atol=1e-13)
            signature=shared_clock_signature(maps,loads)
            self.assertAlmostEqual(signature['inferred_shared_survival'],s,places=12)
            self.assertLess(signature['signature_spread'],1e-12)

    def test_shared_clock_arbitrary_rank_witness(self):
        for rank in range(2,8):
            q=np.zeros(2*rank-1);q[rank-1]=1.;gates=np.ones(len(q)-1)
            for s in (.5,.8,.95):
                block=kinetic_contacts(q,gates,gates,s,True)[:rank,rank:][::-1]
                expected=s**np.maximum.outer(np.arange(rank),np.arange(rank))
                np.testing.assert_allclose(block,expected,atol=1e-14)
                sign,logdet=np.linalg.slogdet(block)
                self.assertEqual(sign,1)
                self.assertAlmostEqual(logdet,rank*(rank-1)/2*np.log(s)+(rank-1)*np.log1p(-s),places=10)


if __name__=='__main__':unittest.main()
