"""
Tests for ProteinGym DMS Benchmark Evaluator
"""

import unittest
from enzyme_gen.proteingym import ProteinGymBenchmarkEvaluator

class TestProteinGym(unittest.TestCase):
    def test_canonical_assay_run(self):
        evaluator = ProteinGymBenchmarkEvaluator()
        res = evaluator.run_assay_evaluation("BLAT_ECOLX_Ostermeier2014", n_mutants=30)
        self.assertEqual(res.target_protein, "TEM-1 Beta-Lactamase")
        self.assertEqual(res.number_of_mutants, 30)
        # Spearman rho should be positive and statistically meaningful
        self.assertGreater(res.spearman_rho, 0.40)
        self.assertEqual(res.benchmark_status, "PASS")

if __name__ == "__main__":
    unittest.main()
