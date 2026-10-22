import unittest
import numpy as np
from clock_identification import identify_clock
from study_kinetic_transfer import kinetic_contacts
from study_load_probes import patterns


class ClockIdentificationTests(unittest.TestCase):
    def test_bounded_band_errors_enclose_clock_and_raw_gates(self):
        rng=np.random.default_rng(81);m=12;loads=patterns(m,'alternating',None)
        for eta in (0.,.0001,.005):
            for _ in range(10):
                gates=rng.uniform(.2,.8,(2,m-1));s=float(rng.uniform(.4,.95))
                maps=[kinetic_contacts(q,*gates,s,True) for q in loads]
                second=np.stack([np.diag(d,2) for d in maps]);third=np.stack([np.diag(d,3) for d in maps])
                fit=identify_clock(loads,second+rng.uniform(-eta,eta,second.shape),third+rng.uniform(-eta,eta,third.shape),eta)
                self.assertFalse(fit['empty']);lo,hi=fit['survival_interval']
                self.assertLessEqual(lo-1e-12,s);self.assertLessEqual(s,hi+1e-12)
                self.assertTrue(np.all(fit['gate_low']-1e-12<=gates));self.assertTrue(np.all(gates<=fit['gate_high']+1e-12))

    def test_non_geometric_family_can_be_rejected(self):
        rng=np.random.default_rng(82);m=16;gates=rng.uniform(.2,.6,(2,m-1));loads=patterns(m,'alternating',None)
        maps=[.5*(kinetic_contacts(q,*gates,.3,True)+kinetic_contacts(q,*gates,.95,True)) for q in loads]
        fit=identify_clock(loads,np.stack([np.diag(d,2) for d in maps]),np.stack([np.diag(d,3) for d in maps]),0.)
        self.assertTrue(fit['empty'])


if __name__=='__main__':unittest.main()
