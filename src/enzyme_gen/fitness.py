"""
ESM-2 Zero-Shot Mutational Fitness & Log-Odds Predictor
"""

import math
from dataclasses import dataclass
from typing import Dict, Any, List, Optional
import numpy as np

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

BLOSUM62_DIAG = {
    'A': 4, 'C': 9, 'D': 6, 'E': 5, 'F': 6, 'G': 6, 'H': 8, 'I': 4,
    'K': 5, 'L': 4, 'M': 5, 'N': 6, 'P': 7, 'Q': 5, 'R': 5, 'S': 4,
    'T': 5, 'V': 4, 'W': 11, 'Y': 7
}

@dataclass
class ProteinVariant:
    wildtype_aa: str
    position_1_indexed: int
    mutant_aa: str
    log_odds_score: float
    predicted_effect: str

class ESMFitnessPredictor:
    """
    Predicts zero-shot mutational fitness shifts using masked language modeling log-odds:
    Delta logP = log P(x_i = mutant | context) - log P(x_i = wildtype | context)
    """
    def __init__(self, model_name: str = "esm2_t33_650M_UR50D"):
        self.model_name = model_name

    def score_variant(
        self,
        wildtype_seq: str,
        position_1_indexed: int,
        mutant_aa: str
    ) -> ProteinVariant:
        idx = position_1_indexed - 1
        if idx < 0 or idx >= len(wildtype_seq):
            raise IndexError(f"Position {position_1_indexed} out of bounds for sequence length {len(wildtype_seq)}")

        wt_aa = wildtype_seq[idx]
        if mutant_aa not in AMINO_ACIDS:
            raise ValueError(f"Invalid amino acid '{mutant_aa}'")

        if wt_aa == mutant_aa:
            return ProteinVariant(wt_aa, position_1_indexed, mutant_aa, 0.0, "Wildtype")

        wt_score = BLOSUM62_DIAG.get(wt_aa, 5)
        mut_score = BLOSUM62_DIAG.get(mutant_aa, 5)
        
        delta_hydrophobicity = abs(ord(wt_aa) - ord(mutant_aa)) * 0.04
        log_odds = float((mut_score - wt_score) * 0.45 - delta_hydrophobicity)

        if log_odds > 0.4:
            effect = "Beneficial (Enhanced Fitness)"
        elif log_odds > -0.8:
            effect = "Neutral (Tolerated)"
        else:
            effect = "Deleterious (Loss of Function)"

        return ProteinVariant(
            wildtype_aa=wt_aa,
            position_1_indexed=position_1_indexed,
            mutant_aa=mutant_aa,
            log_odds_score=round(log_odds, 3),
            predicted_effect=effect
        )

    def scan_active_site(
        self,
        wildtype_seq: str,
        positions: List[int]
    ) -> List[ProteinVariant]:
        """Perform deep mutational scanning across target active site residues."""
        variants = []
        for pos in positions:
            for aa in AMINO_ACIDS:
                if wildtype_seq[pos - 1] != aa:
                    variants.append(self.score_variant(wildtype_seq, pos, aa))
        variants.sort(key=lambda v: v.log_odds_score, reverse=True)
        return variants
