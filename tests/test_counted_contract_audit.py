import json
import unittest
from pathlib import Path

import audit_counted_contract as audit

ROOT = Path(__file__).resolve().parents[1]


class CountedContractAuditTests(unittest.TestCase):
    def test_counterexample_and_proof_checks_pass(self):
        self.assertEqual(audit.audit_counterexample()['gap'], '9/286')
        result = audit.audit_sharpness(trials=500)
        self.assertEqual(result['feasible_points_reconstructed'], 500)
        self.assertTrue(audit.audit_lp()['closed_form_equals_vertex_enumeration'])

    def test_scope_boundaries(self):
        self.assertEqual(audit.audit_false_positive()['kappa'], '98/125')
        self.assertEqual(audit.audit_kappa_above_one()['condition_used_by_proof'], 'kU*ell_max <= 1')
        self.assertIn('survives', audit.audit_state_resolved()['reversible_unknown_generator'])

    def test_saved_report_matches(self):
        saved = json.loads((ROOT/'results/counted_contract_audit.json').read_text())
        self.assertEqual(saved['counterexample'], audit.audit_counterexample())
        self.assertEqual(saved['combined_tolerance_pp'], audit.combined_tolerance())
        self.assertEqual(saved['structural_width_floor_pp'], audit.structural_floor())


if __name__ == '__main__':
    unittest.main()
