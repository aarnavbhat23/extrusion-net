import unittest
import numpy as np
from capture_survival import absence_bounds
from occupancy import capture_snapshot


class CaptureSurvivalTests(unittest.TestCase):
    def test_laplace_identity_matches_independent_full_state_flow(self):
        rng=np.random.default_rng(627)
        for _ in range(160):
            point=np.r_[rng.uniform(.03,1.5),rng.uniform(size=2),rng.uniform(.15,1.2)]
            gamma=float(rng.choice([.1,.3,1.,3.]));lo,hi=absence_bounds(point,point,gamma,2.)
            exact=capture_snapshot(point,gamma,2.)[-1]
            self.assertLessEqual(lo[0],exact+1e-14);self.assertGreaterEqual(hi[0],exact-1e-14)
            self.assertLess(hi[0]-lo[0],1e-12)

    def test_order_corners_include_boxes_spanning_protection_and_harm(self):
        rng=np.random.default_rng(721)
        for _ in range(30):
            low=np.array([.05,.1,.2,rng.uniform(.15,.95)]);high=np.array([1.4,.9,.95,rng.uniform(1.01,1.2)])
            lo,hi=absence_bounds(low,high)
            for _ in range(20):
                point=low+rng.uniform(size=4)*(high-low);value=capture_snapshot(point)[-1]
                self.assertLessEqual(lo[0],value+1e-14);self.assertGreaterEqual(hi[0],value-1e-14)


if __name__=='__main__':unittest.main()
