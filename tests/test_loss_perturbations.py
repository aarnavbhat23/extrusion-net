import unittest
import numpy as np
from study_loss_perturbations import snapshot,recover,DESIGNS


class CompetingLossTests(unittest.TestCase):
    def test_exact_inverse_with_one_sided_capture(self):
        for p in ((.8,1.2,.4,0.,.7,.5,.2),(1.3,.7,.6,.8,0.,1.1,.1)):
            for name in ('three_conditions','absolute_two_turnovers'):
                conditions,conditional,_=DESIGNS[name]
                q=[snapshot(p,c,conditional,True) for c in conditions]
                np.testing.assert_allclose(recover(q,not conditional),p,atol=1e-11,rtol=0)

    def test_no_capture_parameter_is_irrelevant(self):
        p=(1.1,.9,.4,0.,0.,.3,.15)
        conditions,conditional,_=DESIGNS['three_conditions']
        recovered=recover([snapshot(p,c,conditional,True) for c in conditions])
        self.assertIsNone(recovered[5])
        np.testing.assert_allclose(recovered[:5]+recovered[6:],list(p[:5])+list(p[6:]),atol=1e-11,rtol=0)
        other=(*p[:5],1.4,p[6])
        np.testing.assert_allclose(snapshot(p,(1.,3.,.2),False,True),snapshot(other,(1.,3.,.2),False,True),atol=1e-13,rtol=0)

    def test_two_clocks_do_not_resolve_loss_decomposition(self):
        p=(1.1,.9,.4,.4,.6,.5,.15)
        other=(1.1,.9,.35,.4,.6,.15/.35,.2)
        for gamma in (.05,.2,.8):
            np.testing.assert_allclose(snapshot(p,(1.,1.,gamma),False,True),snapshot(other,(1.,1.,gamma),False,True),atol=1e-13,rtol=0)
        self.assertGreater(np.abs(snapshot(p,(1.,2.,.3),False,True)-snapshot(other,(1.,2.,.3),False,True)).sum(),.01)


if __name__ == '__main__': unittest.main()
