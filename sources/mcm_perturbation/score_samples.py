"""Calculate MCM sample program scores using the cached GEO and orthology references."""

import argparse
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from .inputs import read_human_tpm, read_mouse_counts, read_mouse_fpkm, sample_manifest
from .paths import DATA_DIR, FIGURES_DIR, PROJECT_ROOT, REFERENCE_DIR, RESULTS_DIR
from .plotting import save_program_score_figures
from .scoring import calculate_scores, human_symbol_mapping, load_programs, log2_normalized


def tmm_logcpm(counts, manifest, rscript):
    with TemporaryDirectory(prefix="mcm_score_") as work:
        work = Path(work)
        counts.to_csv(work / "counts.tsv", sep="\t")
        manifest[["sample_id", "condition"]].to_csv(work / "samples.tsv", sep="\t", index=False)
        subprocess.run([
            rscript, str(Path(__file__).with_name("count_logcpm.R")),
            str(work / "counts.tsv"), str(work / "samples.tsv"), str(work / "logcpm.tsv"),
        ], check=True)
        result = pd.read_csv(work / "logcpm.tsv", sep="\t", index_col=0)
    if not result.columns.equals(counts.columns) or not result.index.isin(counts.index).all():
        raise ValueError("edgeR output identifiers do not align with the raw counts")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable with edgeR available")
    parser.add_argument("--output", type=Path, default=RESULTS_DIR, help="Directory for result tables")
    parser.add_argument("--figures", type=Path, default=FIGURES_DIR, help="Directory for numbered figures")
    args = parser.parse_args()
    programs = load_programs(PROJECT_ROOT)
    geo = pd.read_csv(REFERENCE_DIR / "geo_samples.tsv", sep="\t", dtype=str)
    mouse_mapping = pd.read_csv(REFERENCE_DIR / "ortholog_symbols.tsv.gz", sep="\t").rename(
        columns={"mouse_symbol": "target_key"}
    )
    score_tables, manifest_tables, coverage_records, mapping_tables = [], [], [], []
    for dataset, recipient in (
        ("GSE175485", "astrocyte"), ("GSE228927", "astrocyte"),
        ("GSE282665", "astrocyte"), ("GSE282665", "motor neuron"),
    ):
        if dataset == "GSE175485":
            raw = read_mouse_counts(DATA_DIR)
        elif dataset == "GSE228927":
            raw = read_mouse_fpkm(DATA_DIR)
        else:
            raw, symbols = read_human_tpm(DATA_DIR, recipient)
        manifest = sample_manifest(geo, dataset, recipient, raw.columns)
        planned = manifest.index[manifest["planned_sample"]]
        if dataset == "GSE175485":
            expression = tmm_logcpm(raw.loc[:, planned], manifest.loc[planned], args.rscript)
            scale = "TMM log2 CPM"
        else:
            expression = log2_normalized(raw, planned)
            scale = "log2(FPKM+1)" if dataset == "GSE228927" else "log2(TPM+1)"
        mapping = human_symbol_mapping(symbols) if dataset == "GSE282665" else mouse_mapping
        manifest["scored"] = manifest.index.isin(expression.columns)
        manifest_tables.append(manifest)
        endpoints = ("SNAP-n",) if recipient == "motor neuron" else ("ACHES", "SNAP-a")
        for program in endpoints:
            values, audit, coverage = calculate_scores(expression, raw.index, programs[program], mapping)
            metadata = {"dataset": dataset, "recipient": recipient, "program": program, "scale": scale}
            score_tables.append(
                manifest.loc[expression.columns].join(values).assign(program=program, scale=scale,
                                                                      low_coverage=coverage["low_coverage"])
            )
            mapping_tables.append(audit.assign(**metadata))
            coverage_records.append({**metadata, **coverage})

    scores = pd.concat(score_tables, ignore_index=True)
    planned_scores = scores.loc[scores["planned_sample"]].copy()
    sample_cols = ["dataset", "recipient", "condition", "sample_id", "gsm_id"]
    wide = planned_scores.pivot(index=sample_cols, columns="program", values="score").reset_index()
    pairs = planned_scores.loc[planned_scores["pair_id"].notna()].pivot(
        index=["dataset", "recipient", "pair_id", "program"], columns="condition", values="score",
    ).reset_index()
    pairs["LPS_minus_resting"] = pairs["LPS-MCM"] - pairs["Resting-MCM"]
    summary = planned_scores.groupby(
        ["dataset", "recipient", "program", "condition"], sort=False,
    )["score"].agg(n_sample_columns="count", mean="mean", sd="std").reset_index()
    if scores.duplicated(subset=["dataset", "recipient", "sample_id", "program"]).any():
        raise ValueError("Duplicate sample program scores")
    args.output.mkdir(parents=True, exist_ok=True)
    outputs = {
        "sample_manifest.tsv": pd.concat(manifest_tables, ignore_index=True),
        "program_scores.tsv": scores, "planned_sample_scores.tsv": wide,
        "program_coverage.tsv": pd.DataFrame(coverage_records), "condition_summary.tsv": summary,
        "paired_score_differences.tsv": pairs,
        "program_mapping.tsv.gz": pd.concat(mapping_tables, ignore_index=True),
    }
    for name, frame in outputs.items():
        frame.to_csv(args.output / name, sep="\t", index=False, float_format="%.12g")
    save_program_score_figures(scores, outputs["program_coverage.tsv"], args.output, args.figures)
    print(f"{len(scores)} score rows; {len(wide)} planned samples -> {args.output}")


if __name__ == "__main__":
    main()
