"""Estimate score differences and pointwise t confidence intervals."""

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, t


CONDITIONS = ("Resting-MCM", "LPS-MCM", "Hemin-MCM")
PROGRAMS = {
    "GSE228927": ("ACHES", "SNAP-a"),
    "GSE282665": ("ACHES", "SNAP-a", "SNAP-n"),
}


def estimate_contrast(scores, treatment, *, paired):
    """Use treatment minus resting, with paired or Welch 95% intervals."""
    reference = CONDITIONS[0]
    rows = scores.loc[scores["condition"].isin([reference, treatment])]
    if rows["sample_id"].duplicated().any() or not np.isfinite(rows["score"]).all():
        raise ValueError("Contrast requires unique sample IDs and finite scores")
    treated = rows.loc[rows["condition"].eq(treatment), "score"].to_numpy()
    control = rows.loc[rows["condition"].eq(reference), "score"].to_numpy()
    if min(len(treated), len(control)) < 2:
        raise ValueError("At least two observations per condition are required")
    if paired:
        if rows["pair_id"].isna().any():
            raise ValueError("A paired contrast requires a pair ID for every sample")
        wide = rows.pivot(index="pair_id", columns="condition", values="score")
        if wide.isna().any().any():
            raise ValueError("Each pair must contain both conditions")
        differences = wide[treatment] - wide[reference]
        effect = differences.mean()
        df = len(differences) - 1
        se = differences.std(ddof=1) / np.sqrt(len(differences))
        method = "paired t"
        n_pairs = len(differences)
    else:
        effect = treated.mean() - control.mean()
        treatment_var = treated.var(ddof=1) / len(treated)
        reference_var = control.var(ddof=1) / len(control)
        se = np.sqrt(treatment_var + reference_var)
        denom = treatment_var**2 / (len(treated) - 1) + reference_var**2 / (len(control) - 1)
        df = (treatment_var + reference_var)**2 / denom if denom > 0 else np.nan
        method = "Welch t"
        n_pairs = None
    estimable = se > 0 and np.isfinite(df)
    half_width = t.ppf(0.975, df) * se if estimable else np.nan
    p_value = 2 * t.sf(abs(effect / se), df) if estimable else np.nan
    return {
        "reference": reference, "treatment": treatment, "method": method,
        "n_reference": len(control), "n_treatment": len(treated), "n_pairs": n_pairs,
        "effect": effect, "standard_error": se, "df": df,
        "ci_low": effect - half_width, "ci_high": effect + half_width, "p_value": p_value,
    }


def summarize_score_effects(scores):
    """Estimate the seven contrasts plotted for the selected datasets."""
    planned = scores.loc[scores["planned_sample"]]
    records = []
    for dataset, programs in PROGRAMS.items():
        paired = dataset == "GSE282665"
        treatments = CONDITIONS[1:2] if paired else CONDITIONS[1:]
        assumption = (
            "Recipient lines matched by GEO metadata; assumes independent line differences."
            if paired else "Exploratory; assumes independent groups. Experiment pairing unavailable."
        )
        for program in programs:
            subset = planned.loc[planned["dataset"].eq(dataset) & planned["program"].eq(program)]
            recipient = subset["recipient"].unique().tolist()
            if len(recipient) != 1:
                raise ValueError(f"Expected one recipient for {dataset} {program}")
            for treatment in treatments:
                records.append({
                    "dataset": dataset, "recipient": recipient[0], "program": program,
                    **estimate_contrast(subset, treatment, paired=paired),
                    "confidence_level": 0.95, "assumption": assumption,
                    "low_coverage": bool(subset["low_coverage"].any()),
                })
    effects = pd.DataFrame(records)
    is_primary = effects["treatment"].eq("LPS-MCM") & effects["recipient"].eq("astrocyte")
    effects["fdr_family"] = np.where(is_primary, "primary", "secondary")
    effects["q_value"] = np.nan
    for _, group in effects.groupby("fdr_family"):
        # Untestable planned comparisons still count toward the family size.
        p = group["p_value"]
        adjusted = false_discovery_control(p.fillna(1).to_numpy(), method="bh")
        effects.loc[group.index, "q_value"] = np.where(p.notna(), adjusted, np.nan)
    return effects
