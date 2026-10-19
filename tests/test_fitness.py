"""
Tests for ESM-2 Zero-Shot Mutational Fitness
"""

import unittest
from enzyme_gen.fitness import ESMFitnessPredictor

class TestFitness(unittest.TestCase):
    def test_single_variant_score(self):
        predictor = ESMFitnessPredictor()
        seq = "MSIQHFRVALIPFFA"
        # Mutate 5th amino acid (H) to A
        v = predictor.score_variant(seq, 5, "A")
        self.assertEqual(v.wildtype_aa, "H")
        self.assertEqual(v.position_1_indexed, 5)
        self.assertEqual(v.mutant_aa, "A")
        self.assertIsInstance(v.log_odds_score, float)
        self.assertIn(v.predicted_effect, [
            "Beneficial (Enhanced Fitness)",
            "Neutral (Tolerated)",
            "Deleterious (Loss of Function)"
        ])

    def test_wildtype_identity(self):
        predictor = ESMFitnessPredictor()
        seq = "MSIQHFRVALIPFFA"
        v = predictor.score_variant(seq, 1, "M")
        self.assertEqual(v.log_odds_score, 0.0)
        self.assertEqual(v.predicted_effect, "Wildtype")

    def test_active_site_scanning(self):
        predictor = ESMFitnessPredictor()
        seq = "MSIQHFRVALIPFFA"
        variants = predictor.scan_active_site(seq, [2, 4])
        # 2 positions * 19 non-wt amino acids = 38 variants
        self.assertEqual(len(variants), 38)
        # Should be sorted descending by log_odds_score
        self.assertGreaterEqual(variants[0].log_odds_score, variants[-1].log_odds_score)

if __name__ == "__main__":
    unittest.main()
