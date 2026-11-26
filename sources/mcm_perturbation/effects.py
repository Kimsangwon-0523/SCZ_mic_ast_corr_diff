"""Estimate score differences and pointwise t confidence intervals."""

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, t


HUMAN_DATASET = "GSE282665"
CONDITIONS = ("Resting-MCM", "LPS-MCM")
PROGRAMS = ("ACHES", "SNAP-a", "SNAP-n")


def estimate_contrast(scores):
    """Estimate LPS minus resting across matched human recipient lines."""
    reference, treatment = CONDITIONS
    rows = scores.loc[scores["condition"].isin(CONDITIONS)]
    if rows["sample_id"].duplicated().any() or not np.isfinite(rows["score"]).all():
        raise ValueError("Contrast requires unique sample IDs and finite scores")
    if rows["pair_id"].isna().any():
        raise ValueError("A paired contrast requires a pair ID for every sample")
    wide = rows.pivot(index="pair_id", columns="condition", values="score")
    if set(wide.columns) != set(CONDITIONS) or wide.isna().any().any() or len(wide) < 2:
        raise ValueError("At least two complete pairs are required")
    differences = wide[treatment] - wide[reference]
    n_pairs = len(differences)
    effect = differences.mean()
    df = n_pairs - 1
    se = differences.std(ddof=1) / np.sqrt(n_pairs)
    estimable = se > 0 and np.isfinite(df)
    half_width = t.ppf(0.975, df) * se if estimable else np.nan
    p_value = 2 * t.sf(abs(effect / se), df) if estimable else np.nan
    return {
        "reference": reference, "treatment": treatment, "method": "paired t",
        "n_reference": n_pairs, "n_treatment": n_pairs, "n_pairs": n_pairs,
        "effect": effect, "standard_error": se, "df": df,
        "ci_low": effect - half_width, "ci_high": effect + half_width, "p_value": p_value,
    }


def summarize_score_effects(scores):
    """Estimate the three GSE282665 program contrasts."""
    planned = scores.loc[scores["planned_sample"] & scores["dataset"].eq(HUMAN_DATASET)]
    records = []
    for program in PROGRAMS:
        subset = planned.loc[planned["program"].eq(program)]
        recipient = "motor neuron" if program == "SNAP-n" else "astrocyte"
        if set(subset["recipient"]) != {recipient}:
            raise ValueError(f"Expected {recipient} scores for {HUMAN_DATASET} {program}")
        records.append({
            "dataset": HUMAN_DATASET, "recipient": recipient, "program": program,
            **estimate_contrast(subset),
            "confidence_level": 0.95,
            "assumption": "Recipient lines matched by GEO metadata; assumes independent line differences.",
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
