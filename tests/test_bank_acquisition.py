import unittest
import numpy as np
from bank_acquisition import confidence_intervals,choose_action


class BankAcquisitionTests(unittest.TestCase):
    def test_multinomial_marginal_intervals_keep_boundary_probabilities(self):
        counts=np.array([0,1000,0]);lo,hi=confidence_intervals(counts)
        np.testing.assert_equal(lo[[0,2]],0.);self.assertEqual(hi[1],1.)
        self.assertTrue(np.all(lo<=hi))

    def test_action_uses_target_distinction_and_no_hidden_observations(self):
        predictions=np.full((3,4,2),.5)
        predictions[1,2]=[.2,.8];predictions[2,2]=[.8,.2]
        targets=np.array([[.1,.9],[.5,.5],[.9,.1]])
        alive=np.ones(3,bool)
        self.assertEqual(choose_action(predictions,targets,alive,'target_disagreement',0),2)
        self.assertEqual(choose_action(predictions,targets,alive,'round_robin',5),1)


if __name__=='__main__':unittest.main()
