import unittest
from fractions import Fraction as F
from counted_loading import counted_bounds, reference_ratio


class CountedLoadingTests(unittest.TestCase):
    def test_identical_reference_and_retention_laws_different_targets(self):
        observed, targets = [], []
        gamma, ell = F(3, 10), F(3, 4)
        for mu, delta, kappa, q1 in [(F(1, 5), F(1, 5), F(4, 5), F(9, 10)),
                                     (F(3, 10), F(1, 10), F(8, 9), F(81, 100))]:
            observed.append([ell*gamma/(gamma+mu+delta), kappa*ell*gamma/(gamma+2*mu+delta),
                             F(9, 10)*ell, q1*kappa*ell])
            targets.append(1-ell*gamma/(gamma+3*mu+delta))
        self.assertEqual(observed[0], observed[1])
        self.assertEqual(observed[0], [F(9, 28), F(1, 5), F(27, 40), F(27, 50)])
        self.assertEqual(targets, [F(35, 44), F(43, 52)])
        self.assertEqual(targets[1]-targets[0], F(9, 286))

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
