"""Repository locations shared by the perturbation commands."""

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data" / "perturbation_validation"
REFERENCE_DIR = DATA_DIR / "reference"
RESULTS_DIR = DATA_DIR / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"
