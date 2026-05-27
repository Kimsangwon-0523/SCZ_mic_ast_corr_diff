import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.stats import norm, gaussian_kde, mannwhitneyu
from typing import Dict

def get_module_eigen_expression(adata, module_dict: Dict[str, list]):
    """
    Calculate average expression (eigengene-like) for each module.
    """
    eigen_df = pd.DataFrame(index=adata.obs_names)
    
    # Normalize var_names for matching
    norm_var_names = pd.Series(adata.var_names).str.replace('.', '-', regex=False).str.strip().tolist()
    var_name_mapping = dict(zip(norm_var_names, adata.var_names))
    
    for mod_name, genes in module_dict.items():
        # Match normalized genes to original var_names
        existing_genes = [var_name_mapping[g] for g in genes if g in var_name_mapping]
        if len(existing_genes) > 0:
            eigen_df[mod_name] = adata[:, existing_genes].X.mean(axis=1)
        else:
            eigen_df[mod_name] = 0
    return eigen_df

def get_cross_corr(eigen1: pd.DataFrame, eigen2: pd.DataFrame):
    """
    Calculate cross-correlation matrix between two sets of module expressions.
    """
    # Pearson correlation
    corr_matrix = pd.DataFrame(index=eigen1.columns, columns=eigen2.columns)
    for col1 in eigen1.columns:
        for col2 in eigen2.columns:
            corr_matrix.loc[col1, col2] = np.corrcoef(eigen1[col1], eigen2[col2])[0, 1]
    return corr_matrix.astype(float)

def fisher_z(r):
    """
    Fisher Z-transformation of correlation coefficient.
    """
    r = np.clip(r, -0.999999, 0.999999)
    return 0.5 * np.log((1 + r) / (1 - r))

def calculate_delta_corr_significance(r1_matrix, n1, r2_matrix, n2):
    """
    Calculate difference in correlations (Fisher Z transformed) and associated p-values.
    """
    z1 = fisher_z(r1_matrix)
    z2 = fisher_z(r2_matrix)
    
    delta_z = z1 - z2
    
    # Standard error of difference between two independent Fisher Z scores
    se = np.sqrt(1/(n1-3) + 1/(n2-3))
    
    z_scores = delta_z / se
    p_values = 2 * (1 - norm.cdf(np.abs(z_scores)))
    
    return delta_z, pd.DataFrame(p_values, index=r1_matrix.index, columns=r1_matrix.columns)

def plot_delta_corr_heatmap(delta_z_df, p_val_df, title="Correlation Difference", 
                           xlabel="Astrocyte modules", ylabel="Microglia modules",
                           vmax=None, cmap="RdBu_r", save_path=None):
    """
    Plot heatmap of correlation differences with significance asterisks.
    """
    n_rows, n_cols = delta_z_df.shape
    
    if vmax is None:
        vmax = np.nanmax(np.abs(delta_z_df.values.flatten()))
    
    norm_scale = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    
    # Create annotation matrix
    annot_matrix = np.full((n_rows, n_cols), "", dtype=object)
    for i in range(n_rows):
        for j in range(n_cols):
            p = p_val_df.iloc[i, j]
            if p < 0.001:
                annot_matrix[i, j] = "***"
            elif p < 0.01:
                annot_matrix[i, j] = "**"
            elif p < 0.05:
                annot_matrix[i, j] = "*"
                
    fig, ax = plt.subplots(figsize=(10, 10))

    delta_z_df = delta_z_df.rename(columns={"Ast15":"Ast15 \n (ACHES)"})
    
    sns.heatmap(
        delta_z_df,
        ax=ax,
        cmap=cmap,
        center=0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.5,
        linecolor="white",
        square=True,
        annot=False,
        cbar_kws={
            "shrink": 0.5,
            "aspect": 20,
            "label": "ΔZ (SCZ - Control)",
            "pad": 0.02,
        },
    )
    
    # Add asterisks
    for i in range(n_rows):
        for j in range(n_cols):
            if annot_matrix[i, j]:
                ax.text(
                    j + 0.92, i + 0.1,
                    annot_matrix[i, j],
                    ha="right", va="top",
                    fontsize=9, fontweight="bold",
                    color="white" if np.abs(delta_z_df.iloc[i, j]) > vmax/2 else "black",
                )

    
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=7, width=0.5, length=3)
    cbar.outline.set_linewidth(0.5)
    cbar.set_label(r"$\Delta Z$  ($Z_{\mathrm{SCZ}} - Z_{\mathrm{Ctrl}}$)", fontsize=10, labelpad=6)
                
    ax.set_title(title, fontsize=16, pad=20)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path+".pdf", dpi=300, bbox_inches="tight")
        plt.savefig(save_path+".png", dpi=600, bbox_inches="tight")
        plt.savefig(save_path+".svg", dpi=600, bbox_inches="tight")

    plt.show()

def plot_raincloud_dz(delta_z_df, target_ast="Ast15", save_path=None):
    """
    Reproduce Figure 2b: Raincloud plot showing ΔZ distribution for 
    pairs involving target_ast vs other pairs, with Mann-Whitney U test.
    """
    # 1. Flatten the matrix and categorize
    df_dz = delta_z_df.copy()
    records = []
    for mic in df_dz.index:
        for ast in df_dz.columns:
            records.append({
                "mic": mic,
                "ast": ast,
                "dz": df_dz.loc[mic, ast],
            })

    flat = pd.DataFrame(records)
    target_dz = flat[flat["ast"] == target_ast]["dz"].values
    other_dz = flat[flat["ast"] != target_ast]["dz"].values

    # 2. Mann-Whitney U test
    stat, p_val = mannwhitneyu(target_dz, other_dz, alternative="greater")
    
    if p_val < 0.001:
        p_str = f"p = {p_val:.2e}"
    else:
        p_str = f"p = {p_val:.4f}"
    
    print(f"Mann-Whitney U = {stat:.1f}, {p_str} (one-sided, greater)")

    # 3. Setup KDE
    x_min = min(target_dz.min(), other_dz.min()) - 0.15
    x_max = max(target_dz.max(), other_dz.max()) + 0.15
    x_grid = np.linspace(x_min, x_max, 400)

    kde_target = gaussian_kde(target_dz, bw_method=0.4)(x_grid)
    kde_other = gaussian_kde(other_dz, bw_method=0.3)(x_grid)

    kde_target_norm = kde_target / kde_target.max() * 0.35
    kde_other_norm = kde_other / kde_other.max() * 0.35

    # 4. Figure Layout
    plt.rcParams.update({
        "font.family": "Liberation Serif",
        "font.size": 9,
        "mathtext.default": "regular",
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
    })

    fig, ax = plt.subplots(figsize=(7, 4))
    np.random.seed(42)
    
    y_other = 0
    y_target = 1

    # --- Other pairs (bottom) ---
    ax.fill_between(
        x_grid, y_other, y_other + kde_other_norm,
        color="#D4D4D4", alpha=0.7, zorder=2,
    )
    ax.plot(x_grid, y_other + kde_other_norm, color="#A0A0A0", linewidth=0.8, zorder=2)

    ax.boxplot(
        other_dz, positions=[y_other], vert=False,
        widths=0.08, patch_artist=True,
        boxprops=dict(facecolor="white", edgecolor="#888888", linewidth=0.7),
        medianprops=dict(color="#555555", linewidth=1.0),
        whiskerprops=dict(color="#888888", linewidth=0.7),
        capprops=dict(color="#888888", linewidth=0.7),
        flierprops=dict(marker="", markersize=0),
        zorder=5,
    )

    jitter_other = np.random.uniform(-0.18, -0.05, size=len(other_dz))
    ax.scatter(
        other_dz, y_other + jitter_other,
        s=18, alpha=0.7, color="#B0B0B0", edgecolors="none", zorder=3,
    )

    # --- Target pairs (top) ---
    ax.fill_between(
        x_grid, y_target, y_target + kde_target_norm,
        color="#7B4FB5", alpha=0.7, zorder=2,
    )
    ax.plot(x_grid, y_target + kde_target_norm, color="#5C3A8C", linewidth=0.8, zorder=2)

    ax.boxplot(
        target_dz, positions=[y_target], vert=False,
        widths=0.08, patch_artist=True,
        boxprops=dict(facecolor="white", edgecolor="#5C3A8C", linewidth=0.7),
        medianprops=dict(color="#3D2066", linewidth=1.0),
        whiskerprops=dict(color="#5C3A8C", linewidth=0.7),
        capprops=dict(color="#5C3A8C", linewidth=0.7),
        flierprops=dict(marker="", markersize=0),
        zorder=5,
    )

    jitter_target = np.random.uniform(-0.18, -0.05, size=len(target_dz))
    ax.scatter(
        target_dz, y_target + jitter_target,
        s=18, alpha=0.7, color="#7B4FB5", edgecolors="white", linewidths=0.3, zorder=3,
    )

    # --- Significance bracket ---
    bracket_x = max(target_dz.max(), other_dz.max()) + 0.23
    bracket_y_bottom = y_other
    bracket_y_top = y_target
    bracket_tip = 0.02

    ax.plot([bracket_x, bracket_x], [bracket_y_bottom, bracket_y_top], color="#444444", linewidth=0.8, clip_on=False, zorder=6)
    ax.plot([bracket_x - bracket_tip, bracket_x], [bracket_y_bottom, bracket_y_bottom], color="#444444", linewidth=0.8, clip_on=False, zorder=6)
    ax.plot([bracket_x - bracket_tip, bracket_x], [bracket_y_top, bracket_y_top], color="#444444", linewidth=0.8, clip_on=False, zorder=6)

    ax.text(
        bracket_x + 0.02,
        (bracket_y_bottom + bracket_y_top) / 2,
        f"Mann-Whitney U\n(one-sided)\n{p_str}",
        fontsize=9, color="#333333",
        ha="left", va="center",
        linespacing=1.4,
    )

    # --- Formatting ---
    ax.axvline(0, color="#888888", linewidth=0.5, linestyle="-", alpha=0.3, zorder=0)
    ax.set_yticks([y_other, y_target])
    if target_ast == "Ast15":
        target_ast = "Ast15 (ACHES)"
    ax.set_yticklabels(
        [f"Other module pairs\n(n={len(other_dz)})",
         f"{target_ast} – Mic pairs\n(n={len(target_dz)})"],
        fontsize=11,
    )
    ax.set_xlabel(r"$\Delta Z$  ($Z_{\mathrm{SCZ}} - Z_{\mathrm{Ctrl}}$)", fontsize=12)
    ax.set_ylim(-0.35, 1.55)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)

    plt.tight_layout()
    if save_path:
        for ext in [".pdf", ".png", ".svg"]:
            plt.savefig(save_path + ext, dpi=300, bbox_inches="tight", pad_inches=0.15)
    plt.show()
