from decimal import localcontext, Decimal
import unittest
from study_capture_identifiability import audit, recover
from verify_occupancy_witnesses import decimal_snapshot


class CaptureIdentifiabilityTests(unittest.TestCase):
    def test_inverse_across_capture_and_removal_boundaries(self):
        with localcontext() as context:
            context.prec = 70
            for p in ((.03, 0, 0, .15), (.3, 0, .7, 1.2), (.7, 1, 0, .5), (1.5, 1, 1, 1.2), (.4, 1e-6, 0, .7)):
                for gamma in (.1, .3):
                    audit(p, gamma)

    def test_capture_parameter_irrelevant_only_when_no_capture(self):
        with localcontext() as context:
            context.prec = 70
            p = (.3, .6, .4, .5); changed = (.3, .6, .4, 1.1)
            first, second = decimal_snapshot(p), decimal_snapshot(changed)
            self.assertGreater(abs(first[27]-second[27]), Decimal('.0001'))
            self.assertGreater(abs(recover(first, .3)[-1]-recover(second, .3)[-1]), Decimal('.5'))


if __name__ == '__main__':
    unittest.main()
