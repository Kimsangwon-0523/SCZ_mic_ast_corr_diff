"""Project fixed, signed programs onto sample-level expression."""

import numpy as np
import pandas as pd


def load_programs(repo):
    data = repo / "data"
    aches_path = data / "primary_cohort/WGCNA_modules/modules_astrocyte/astrocyte_module_15.txt"
    aches = aches_path.read_text().split()
    if len(aches) != 363 or len(set(aches)) != 363:
        raise ValueError("ACHES must have 363 distinct source symbols")
    # Resolve the historical dot spelling in this specific source list.
    programs = {"ACHES": pd.Series(1.0, index=pd.Index(aches).map(
        lambda symbol: "HLA-E" if symbol == "HLA.E" else symbol
    ))}
    for program, supplement in (("SNAP-a", 10), ("SNAP-n", 11)):
        path = data / f"SNAP_loadings/41586_2024_7109_MOESM{supplement}_ESM.xlsx"
        weights = pd.read_excel(path, sheet_name="Gene loadings", usecols=["Gene", program])
        weights = weights.set_index("Gene")[program].astype(float)
        if weights.index.has_duplicates or weights.index.isna().any() or not np.isfinite(weights).all():
            raise ValueError(f"Invalid source loadings for {program}")
        programs[program] = weights.loc[weights.ne(0)]
    return programs


def log2_normalized(expression, planned_samples):
    """Apply an all-zero filter based on the planned samples, then log2(x + 1)."""
    values = expression.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Expected finite nonnegative normalized expression")
    return np.log2(expression.loc[expression.loc[:, planned_samples].gt(0).any(axis=1)] + 1)


def human_symbol_mapping(symbols):
    table = symbols.rename("human_symbol").reset_index()
    table = table.loc[table["human_symbol"].ne("")]
    table = table.drop_duplicates("human_symbol", keep=False)
    return table.rename(columns={"gene_id": "target_key"})


def calculate_scores(expression, all_genes, weights, mapping):
    """Return scores and per-gene mapping diagnostics for one positive-mass program."""
    if expression.index.has_duplicates or expression.columns.has_duplicates:
        raise ValueError("Expression must have unique gene keys and sample columns")
    weights = weights.rename("weight").rename_axis("human_symbol").reset_index()
    audit = weights.merge(mapping, on="human_symbol", how="left", validate="one_to_one")
    audit["matrix_present"] = audit["target_key"].isin(all_genes)
    audit["retained"] = audit["target_key"].isin(expression.index)
    audit["status"] = np.select(
        [audit["target_key"].isna(), ~audit["matrix_present"], ~audit["retained"]],
        ["unmapped_or_ambiguous", "absent_from_matrix", "removed_by_expression_filter"],
        default="retained",
    )
    used = audit.loc[audit["retained"]]
    if used.empty or used["target_key"].duplicated().any():
        raise ValueError("Program has no retained genes or conflicting gene mappings")
    mass = used["weight"].abs().sum()
    if mass <= 0:
        raise ValueError("Program has zero retained absolute weight")
    x = expression.loc[used["target_key"]].to_numpy(dtype=float)
    result = pd.DataFrame(
        {"weighted_sum": used["weight"].to_numpy() @ x},
        index=expression.columns.rename("sample_id"),
    )
    result["score"] = result["weighted_sum"] / mass
    coverage = {
        "source_genes": len(audit),
        "mapped_genes": int(audit["target_key"].notna().sum()),
        "matrix_genes": int(audit["matrix_present"].sum()),
        "retained_genes": len(used),
        "gene_coverage": len(used) / len(audit),
        "weight_mass_coverage": mass / audit["weight"].abs().sum(),
    }
    for sign, rows in (("positive", audit["weight"].gt(0)), ("negative", audit["weight"].lt(0))):
        total = audit.loc[rows, "weight"].abs().sum()
        retained_mass = audit.loc[rows & audit["retained"], "weight"].abs().sum()
        coverage[f"{sign}_mass_coverage"] = retained_mass / total if total > 0 else np.nan
    coverage["low_coverage"] = (
        coverage["retained_genes"] < 20
        or min(coverage["gene_coverage"], coverage["weight_mass_coverage"]) < 0.5
    )
    return result, audit, coverage
