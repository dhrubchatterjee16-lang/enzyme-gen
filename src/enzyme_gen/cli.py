"""
Command-Line Interface for EnzymeGen
"""

import argparse
from .fitness import ESMFitnessPredictor
from .proteingym import ProteinGymBenchmarkEvaluator

def main():
    parser = argparse.ArgumentParser(
        prog="enzyme-gen",
        description="EnzymeGen: ProteinGym DMS Benchmark & ESM-2 Zero-Shot Mutational Engine"
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # score subcommand
    score_parser = subparsers.add_parser("score", help="Score single-point variant zero-shot fitness")
    score_parser.add_argument("--sequence", type=str, default="MSIQHFRVALIPFFAAFCLPVFA", help="Wildtype sequence")
    score_parser.add_argument("--pos", type=int, default=5, help="1-indexed position")
    score_parser.add_argument("--mut", type=str, default="A", help="Mutant amino acid")

    # scan subcommand
    scan_parser = subparsers.add_parser("scan", help="Run active site deep mutational scanning")
    scan_parser.add_argument("--sequence", type=str, default="MSIQHFRVALIPFFAAFCLPVFA")
    scan_parser.add_argument("--positions", type=int, nargs="+", default=[3, 5, 7])

    # benchmark subcommand
    bench_parser = subparsers.add_parser("benchmark", help="Run ProteinGym DMS validation assay")
    bench_parser.add_argument("--assay", type=str, default="BLAT_ECOLX_Ostermeier2014")
    bench_parser.add_argument("--mutants", type=int, default=50)

    args = parser.parse_args()

    if args.command == "score":
        predictor = ESMFitnessPredictor()
        v = predictor.score_variant(args.sequence, args.pos, args.mut)
        print("\n--- ESM-2 Zero-Shot Variant Score ---")
        print(f"Variant:         {v.wildtype_aa}{v.position_1_indexed}{v.mutant_aa}")
        print(f"Log-Odds Score:  {v.log_odds_score}")
        print(f"Fitness Effect:  {v.predicted_effect}\n")

    elif args.command == "scan":
        predictor = ESMFitnessPredictor()
        variants = predictor.scan_active_site(args.sequence, args.positions)
        print(f"\n--- In-Silico Deep Mutational Scanning ({len(variants)} variants) ---")
        print(f"{'Rank':<5} {'Variant':<10} {'Log-Odds':<10} {'Effect'}")
        print("-" * 55)
        for i, v in enumerate(variants[:10], 1):
            print(f"{i:<5} {v.wildtype_aa + str(v.position_1_indexed) + v.mutant_aa:<10} {v.log_odds_score:<10} {v.predicted_effect}")
        print()

    elif args.command == "benchmark":
        evaluator = ProteinGymBenchmarkEvaluator()
        res = evaluator.run_assay_evaluation(args.assay, args.mutants)
        print("\n--- ProteinGym DMS Benchmark Result ---")
        print(f"Assay ID:        {res.assay_name}")
        print(f"Target Protein:  {res.target_protein}")
        print(f"Mutants Tested:  {res.number_of_mutants}")
        print(f"Spearman's Rho:  {res.spearman_rho} (p-value: {res.p_value:.2e})")
        print(f"SOTA Baseline:   {res.baseline_sota_rho}")
        print(f"Benchmark Status: {res.benchmark_status}\n")

    else:
        parser.print_help()

if __name__ == "__main__":
    main()
