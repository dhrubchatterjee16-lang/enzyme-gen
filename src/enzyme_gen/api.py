"""
FastAPI Microservice for EnzymeGen
"""

from typing import Dict, Any, List, Optional
from pydantic import BaseModel, Field

try:
    from fastapi import FastAPI, HTTPException
    HAS_FASTAPI = True
except ImportError:
    HAS_FASTAPI = False

from .fitness import ESMFitnessPredictor
from .proteingym import ProteinGymBenchmarkEvaluator

if HAS_FASTAPI:
    app = FastAPI(
        title="EnzymeGen API",
        description="ProteinGym Benchmark & ESM-2 Zero-Shot Mutational Engine",
        version="1.0.0"
    )

    class VariantScoreRequest(BaseModel):
        wildtype_seq: str = Field(..., description="Wildtype amino acid sequence")
        position_1_indexed: int = Field(..., ge=1)
        mutant_aa: str = Field(..., max_length=1)

    class DMSScanRequest(BaseModel):
        wildtype_seq: str
        target_positions: List[int]

    @app.get("/")
    def root():
        return {
            "service": "EnzymeGen",
            "version": "1.0.0",
            "status": "operational"
        }

    @app.post("/variant/score")
    def api_score(req: VariantScoreRequest):
        predictor = ESMFitnessPredictor()
        try:
            v = predictor.score_variant(req.wildtype_seq, req.position_1_indexed, req.mutant_aa)
            return {
                "variant": f"{v.wildtype_aa}{v.position_1_indexed}{v.mutant_aa}",
                "log_odds_score": v.log_odds_score,
                "predicted_effect": v.predicted_effect
            }
        except Exception as e:
            raise HTTPException(status_code=400, detail=str(e))

    @app.post("/variant/scan")
    def api_scan(req: DMSScanRequest):
        predictor = ESMFitnessPredictor()
        variants = predictor.scan_active_site(req.wildtype_seq, req.target_positions)
        return {
            "num_variants": len(variants),
            "top_candidates": [
                {
                    "variant": f"{v.wildtype_aa}{v.position_1_indexed}{v.mutant_aa}",
                    "log_odds_score": v.log_odds_score,
                    "predicted_effect": v.predicted_effect
                }
                for v in variants[:15]
            ]
        }
else:
    app = None
