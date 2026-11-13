import unittest
import numpy as np
from occupancy import capture_snapshot,markov_snapshot,external_snapshots


class OccupancyTests(unittest.TestCase):
    def test_published_forward_and_numpy_resolvent_agree(self):
        from dlem_comparator import load_model,snapshot
        rng=np.random.default_rng(901);model=load_model()
        for _ in range(8):
            parameters=np.r_[rng.uniform(.05,2.,16),rng.uniform(.05,1.)]
            for gamma in (.1,.3,3.):
                np.testing.assert_allclose(snapshot(parameters,gamma,model=model),markov_snapshot(parameters,gamma),atol=1e-13)

    def test_no_capture_reduces_to_markov_and_probability_is_conserved(self):
        for mu in (.15,.35,.7):
            parameters=np.r_[np.ones(16),mu]
            for gamma in (.1,.3,3.):
                np.testing.assert_allclose(capture_snapshot([mu,0.,0.,.25],gamma),markov_snapshot(parameters,gamma),atol=1e-13)
                prediction=capture_snapshot([mu,.8,.5,.25],gamma)
                self.assertTrue(np.all(prediction>=0));self.assertAlmostEqual(prediction.sum(),1.)
                # With constant removal, survival to an independent Exp(gamma) clock is gamma/(gamma+mu).
                self.assertAlmostEqual(capture_snapshot([mu,.8,.5,1.],gamma)[-1],mu/(gamma+mu))

    def test_external_categorical_trials_match_hidden_state_model(self):
        parameters=[.35,.8,.5,.25];prediction=capture_snapshot(parameters)
        observed=external_snapshots(parameters,.3,10000,410000000)/10000
        # Uniform Hoeffding envelope for these 82 categorical events.
        tolerance=np.sqrt(np.log(2*82/1e-5)/(2*10000))
        self.assertLess(np.max(np.abs(observed-prediction)),tolerance)


if __name__=='__main__':unittest.main()
