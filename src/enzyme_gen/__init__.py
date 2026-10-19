"""
EnzymeGen: ProteinGym Deep Mutational Scanning (DMS) Benchmark & ESM-2 Fitness Engine
"""

from .fitness import ESMFitnessPredictor, ProteinVariant
from .proteingym import ProteinGymBenchmarkEvaluator, DMSBenchmarkResult

__version__ = "1.0.0"
__all__ = [
    "ESMFitnessPredictor",
    "ProteinVariant",
    "ProteinGymBenchmarkEvaluator",
    "DMSBenchmarkResult",
    "__version__"
]
