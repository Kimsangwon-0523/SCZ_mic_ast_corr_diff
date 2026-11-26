"""Rebuild the human astrocyte figure and contrast table from saved scores."""

import argparse
from pathlib import Path

import matplotlib
import pandas as pd

from .paths import FIGURES_DIR, RESULTS_DIR

matplotlib.use("Agg")
from .plotting import save_program_score_figures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=RESULTS_DIR, help="Directory for result tables")
    parser.add_argument("--figures", type=Path, default=FIGURES_DIR, help="Directory for numbered figures")
    args = parser.parse_args()
    results = args.results
    scores = pd.read_csv(results / "program_scores.tsv", sep="\t")
    save_program_score_figures(scores, results, args.figures)
    print(f"Saved figure 24 to {args.figures} and program_score_effects_ci.tsv to {results}")


if __name__ == "__main__":
    main()
