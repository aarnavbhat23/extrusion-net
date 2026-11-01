import unittest
from counted_loading import counted_bounds, reference_ratio


class CountedLoadingTests(unittest.TestCase):
    def test_exact_ratio_and_sensitivity_image(self):
        self.assertEqual(reference_ratio([(.6, .6), (.48, .48)]), ((.8, .8), 'ok'))
        ratio, status = reference_ratio([(.6, .6), (.432, .432)], (.9, 1.))
        self.assertEqual(status, 'ok')
        self.assertAlmostEqual(ratio[0], .72)
        self.assertAlmostEqual(ratio[1], .8)

    def test_zero_and_incompatible_references(self):
        self.assertEqual(reference_ratio([(0, .5), (.2, .3)])[1], 'zero_reference_lower')
        self.assertEqual(reference_ratio([(.3, .4), (.5, .6)])[1], 'incompatible_reference')
        result = counted_bounds([30, 20], [100]*2, [0, 0], [100]*2)
        self.assertEqual(result['interval'], [0., 1.])

    def test_uncertainty_cannot_shrink_from_wider_sensitivity_contract(self):
        tight = counted_bounds([7714, 4800], [24000]*2, [18360, 16524], [24000]*2)
        broad = counted_bounds([7714, 4800], [24000]*2, [18360, 16524], [24000]*2, (.9, 1.))
        self.assertEqual(tight['status'], 'ok')
        self.assertEqual(broad['status'], 'ok')
        self.assertLessEqual(broad['interval'][0], tight['interval'][0]+1e-9)
        self.assertGreaterEqual(broad['interval'][1], tight['interval'][1]-1e-9)
