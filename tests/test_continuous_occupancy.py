import unittest
import numpy as np
from continuous_occupancy import interval_snapshot,enclose,DOMAIN_LOW,DOMAIN_HIGH
from occupancy import capture_snapshot


class ContinuousOccupancyTests(unittest.TestCase):
    def test_interval_flow_contains_points_and_refinement_is_nested(self):
        rng=np.random.default_rng(671)
        for _ in range(18):
            endpoints=rng.uniform(size=(2,4))*(DOMAIN_HIGH-DOMAIN_LOW)+DOMAIN_LOW
            lo=endpoints.min(axis=0);hi=endpoints.max(axis=0);mid=(lo+hi)/2
            lower,upper=interval_snapshot(lo,hi,gamma=.1,removal_factor=2.)
            for point in (lo,hi,mid,lo+rng.uniform(size=4)*(hi-lo)):
                truth=capture_snapshot(point,.1,2.)
                self.assertTrue(np.all(lower[0]<=truth+1e-14));self.assertTrue(np.all(upper[0]>=truth-1e-14))
            child_lower,child_upper=interval_snapshot(lo,mid,.1,2.)
            self.assertTrue(np.all(child_lower>=lower-1e-13));self.assertTrue(np.all(child_upper<=upper+1e-13))

    def test_point_boxes_are_precise_and_resource_exhaustion_never_discards_domain(self):
        point=np.array([.35,.2,.8,.25]);lo,hi=interval_snapshot(point,point)
        exact=capture_snapshot(point)
        self.assertTrue(np.all(lo[0]<=exact+1e-14));self.assertTrue(np.all(hi[0]>=exact-1e-14))
        self.assertLess(np.max(hi-lo),1e-12)
        counts=np.rint(exact*100000).astype(int);counts[-1]+=100000-counts.sum()
        result,lower,upper=enclose([(.3,counts)],max_boxes=1)
        self.assertFalse(result['rejected']);self.assertFalse(result['attained'])
        np.testing.assert_equal(lower,[DOMAIN_LOW]);np.testing.assert_equal(upper,[DOMAIN_HIGH])


if __name__=='__main__':unittest.main()
