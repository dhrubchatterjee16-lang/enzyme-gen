# EnzymeGen: ProteinGym DMS Benchmark & ESM-2 Zero-Shot Mutational Engine

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Python: 3.10+](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue.svg)](https://www.python.org/)
[![Hugging Face Spaces](https://img.shields.io/badge/%F0%9F%A4%97%20Hugging%20Face-enzyme-yellow)](https://huggingface.co/spaces/dhrubchatterjee/enzyme)
[![Tests: Passing](https://img.shields.io/badge/tests-4%2F4%20passed-brightgreen.svg)](tests/)

**EnzymeGen** is a computational protein engineering and directed evolution engine. It pairs ESM-2 (Evolutionary Scale Modeling) zero-shot masked language modeling with the gold-standard **ProteinGym** Deep Mutational Scanning (DMS) benchmark suite (217 experimental DMS assays).

---

## 🔬 Scientific Foundations & Mathematical Formulations

### 1. ESM-2 Zero-Shot Log-Odds Fitness Scoring
Zero-shot fitness evaluation evaluates whether a single-point or combinatorial amino acid substitution enhances catalytic stability or causes loss of function:

$$\Delta \log P(x_{\mathrm{mut}} | \mathbf{x}) = \log P\left(x_i = \mathrm{mutant} \mid \mathbf{x}_{-i}\right) - \log P\left(x_i = \mathrm{wildtype} \mid \mathbf{x}_{-i}\right)$$

where $\mathbf{x}_{-i}$ represents the surrounding sequence context masked at residue position $i$.

### 2. ProteinGym DMS Benchmark Validation
Predictions are systematically benchmarked against curated high-throughput experimental mutational landscapes (TEM-1 Beta-Lactamase, GFP, PABP) via Spearman's rank correlation coefficient $\rho$:

$$\rho = 1 - \frac{6 \sum_{i=1}^n d_i^2}{n(n^2 - 1)}$$

where $d_i = \mathrm{rank}(y_i) - \mathrm{rank}(\hat{y}_i)$.

---

## ⚡ Quickstart

### Installation

```bash
pip install -e .
```

### Python API

```python
from enzyme_gen import ESMFitnessPredictor, ProteinGymBenchmarkEvaluator

# 1. Zero-Shot Single-Point Variant Fitness
predictor = ESMFitnessPredictor()
variant = predictor.score_variant(
    wildtype_seq="MSIQHFRVALIPFFAAFCLPVFA",
    position_1_indexed=5,
    mutant_aa="A"
)
print(f"Variant: {variant.wildtype_aa}{variant.position_1_indexed}{variant.mutant_aa}")
print(f"Log-Odds Score: {variant.log_odds_score}")
print(f"Predicted Effect: {variant.predicted_effect}")

# 2. In-Silico Deep Mutational Scanning of Active Site
candidates = predictor.scan_active_site("MSIQHFRVALIPFFAAFCLPVFA", positions=[3, 5, 7])
for c in candidates[:5]:
    print(f"Top Variant: {c.wildtype_aa}{c.position_1_indexed}{c.mutant_aa} -> {c.log_odds_score}")

# 3. ProteinGym Benchmark Validation
evaluator = ProteinGymBenchmarkEvaluator()
res = evaluator.run_assay_evaluation("BLAT_ECOLX_Ostermeier2014", n_mutants=50)
print(f"Assay: {res.target_protein} | Spearman rho: {res.spearman_rho} (Status: {res.benchmark_status})")
```

---

## 🖥️ Command-Line Interface (CLI)

```bash
# Score a specific mutation
enzyme-gen score --sequence "MSIQHFRVALIPFFA" --pos 5 --mut A

# Deep mutational scan across active site residues
enzyme-gen scan --sequence "MSIQHFRVALIPFFA" --positions 3 5 7

# Run ProteinGym benchmark assay
enzyme-gen benchmark --assay BLAT_ECOLX_Ostermeier2014 --mutants 50
```

---

## 🧪 Automated Verification Suite

```bash
python run_tests.py
```
```text
Ran 4 tests in 0.003s
OK - ALL ENZYME-GEN TESTS PASSED (4/4)!
```

---

## 📚 References
- Notin, P., et al. (2023). *ProteinGym: Large-Scale Benchmarks for Protein Fitness Prediction and Design*. Advances in Neural Information Processing Systems (NeurIPS Datasets and Benchmarks).
- Lin, Z., et al. (2023). *Evolutionary-scale prediction of atomic-level protein structure with a language model*. Science, 379(6637), 1123-1130.
- Ostermeier, M., et al. (2014). *Comprehensive mapping of fitness landscapes in TEM-1 beta-lactamase*. Journal of Molecular Biology.

---

## 📄 License
MIT License. Copyright (c) 2026 Dhrub Chatterjee.
