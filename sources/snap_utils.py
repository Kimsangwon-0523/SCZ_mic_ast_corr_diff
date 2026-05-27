import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, t as t_dist
from scipy import sparse
import os

def normalize_gene_name(g: str) -> str:
    """Normalize gene name: '.' -> '-' conversion."""
    return str(g).replace('.', '-').strip()

def get_module_means(pb_ct_df, module_dict, mod_prefix='Ast'):
    """
    Calculate mean expression per module for given pseudobulk data.
    """
    # Normalize column names
    pb_ct_df.columns = [normalize_gene_name(c) for c in pb_ct_df.columns]

    module_means = {}
    for mid, genes in module_dict.items():
        orig_genes = [str(g) for g in genes]
        norm_genes = [normalize_gene_name(g) for g in orig_genes]

        # Unique normalized genes
        norm_genes_unique = list(dict.fromkeys(norm_genes))

        # Match columns
        matched_cols = [g for g in norm_genes_unique if g in pb_ct_df.columns]

        if len(matched_cols) == 0:
            module_means[f"{mod_prefix}{mid}"] = pd.Series(np.nan, index=pb_ct_df.index)
        else:
            module_means[f"{mod_prefix}{mid}"] = pb_ct_df[matched_cols].mean(axis=1, skipna=True)

        print(f"[module {mid}] used {len(matched_cols)}/{len(orig_genes)} genes.")

    return pd.DataFrame(module_means, index=pb_ct_df.index)

def snap_score_calc(adata, loadings):
    """
    Calculate SNAP score by projecting loadings onto expression matrix.
    """
    common_genes = adata.var_names.intersection(loadings.index)
    load_vec = loadings.loc[common_genes].values

    expr = adata[:, common_genes].X
    if sparse.issparse(expr):
        scores = expr.dot(load_vec)
    else:
        scores = np.dot(expr, load_vec)

    adata.obs['SNAP_a_score'] = np.array(scores).ravel()

def regression_ci(x, y, x_pred, ci=0.95):
    """
    Calculate regression line and its 95% confidence interval.
    """
    n = len(x)
    m, b = np.polyfit(x, y, 1)
    y_pred = m * x_pred + b
    y_hat = m * x + b
    resid = y - y_hat
    se = np.sqrt(np.sum(resid**2) / (n - 2))
    x_mean = np.mean(x)
    sx = np.sum((x - x_mean)**2)
    t_val = t_dist.ppf((1 + ci) / 2, df=n - 2)
    ci_band = t_val * se * np.sqrt(1/n + (x_pred - x_mean)**2 / sx)
    return y_pred, y_pred - ci_band, y_pred + ci_band

def plot_snap_ast_scatter(df, x_col, y_col, condition_col, stats, palette, title=None, output_path=None):
    """
    Plot SNAP score vs module score scatter plot with regression and CIs.
    """
    # Style setup
    plt.rcParams.update({
        "font.family": "Liberation Serif",
        "font.size": 9,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
    })

    fig, ax = plt.subplots(figsize=(5, 5))

    conditions = df[condition_col].unique()
    # Sort to maintain consistency (Control first usually)
    conditions = sorted(conditions, key=lambda x: 0 if x == "Control" else 1)

    for cond in conditions:
        sub = df[df[condition_col] == cond]
        c = palette.get(cond, "gray")
        x = sub[x_col].values
        y = sub[y_col].values

        # Scatter
        ax.scatter(
            x, y,
            c=c, s=22, alpha=0.80,
            linewidths=0.2, label=cond, zorder=2,
        )

        # Regression + 95% CI
        x_range = np.linspace(x.min(), x.max(), 200)
        y_fit, ci_lo, ci_hi = regression_ci(x, y, x_range)

        ax.plot(x_range, y_fit, color=c, lw=1, alpha=0.3, zorder=1)
        ax.fill_between(
            x_range, ci_lo, ci_hi,
            color=c, alpha=0.35,
            edgecolor=c, linewidth=0,
            zorder=1,
        )

    # Stat annotation
    stat_texts = []
    for cond in conditions:
        r, p = stats[cond]
        stat_texts.append(f"{cond}: r = {r:.3f}, p = {p:.1e}")
    
    stat_text = "\n".join(stat_texts)
    ax.text(
        0.97, 0.97, stat_text,
        transform=ax.transAxes, ha="right", va="top",
        fontsize=9, 
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.8", alpha=1),
    )

    # Axis labels
    ax.set_xlabel("SNAP-a score", fontsize=10, labelpad=6)
    if y_col == "Ast15":
        y_col = "Ast15 (ACHES)"
    ax.set_ylabel(f"{y_col} score", fontsize=10, labelpad=6)
    if title is None:
        title = f"{y_col} score vs. SNAP-a score by condition"

    ax.set_title(
        title,
        fontsize=14, fontweight="bold", pad=10,
    )

    # Legend
    leg = ax.legend(
        fontsize=9, frameon=True, fancybox=True,
        edgecolor="0.7", framealpha=1,
        loc="lower left", handletextpad=0.4,
        borderpad=0.2,
    )
    if leg:
        leg.get_frame().set_boxstyle("round,pad=0.4")
        leg.get_frame().set_linewidth(0.5)

    # Clean up spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)

    plt.tight_layout()
    if output_path:
        for ext in ['.pdf', '.png', '.svg']:
            plt.savefig(f"{output_path}{ext}", dpi=300, bbox_inches="tight")
    plt.show()

def plot_correlation_shift_arrows(z_ctrl_series, z_scz_series, title, output_path=None):
    """
    Plot arrow plot (lollipop-style shift) showing correlation changes 
    from Control to SCZ across modules.
    """
    import matplotlib.patches as mpatches

    # 1. Prepare DataFrame
    df = pd.DataFrame({
        "module": z_ctrl_series.index,
        "z_ctrl": z_ctrl_series.values,
        "z_scz":  z_scz_series.loc[z_ctrl_series.index].values,
    })
    df["dz"] = df["z_scz"] - df["z_ctrl"]
    df = df.sort_values("z_ctrl", ascending=False).reset_index(drop=True)

    # 2. Colors by direction and intensity
    dz_abs_max = df["dz"].abs().max()
    if dz_abs_max == 0: dz_abs_max = 1
    colors = []
    for v in df["dz"]:
        intensity = 0.3 + 0.6 * (abs(v) / dz_abs_max)
        if v >= 0:
            colors.append(plt.cm.Reds(intensity))
        else:
            colors.append(plt.cm.Blues(intensity))

    # 3. Style setup
    plt.rcParams.update({
        "font.family": "Liberation Serif",
        "font.size": 9,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
    })

    # 4. Plot
    fig, ax = plt.subplots(figsize=(3, 5))
    y_pos = np.arange(len(df))

    for i, row in df.iterrows():
        z0 = row["z_ctrl"]
        z1 = row["z_scz"]
        c  = colors[i]

        # Arrow line (ctrl → scz)
        ax.annotate(
            "",
            xy=(z1, i),          # arrowhead at SCZ
            xytext=(z0, i),      # tail at Control
            arrowprops=dict(
                arrowstyle="-|>",
                color=c,
                lw=1.8,
                mutation_scale=10,
                shrinkA=0,
                shrinkB=0,
            ),
        )

        # Dot at control (start)
        ax.plot(z0, i, "o", color=c, markersize=5, zorder=3,
                markeredgecolor="white", markeredgewidth=0.4)

    # 5. Formatting
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["module"], fontsize=9)
    ax.invert_yaxis()
    ax.axvline(0, color="0.4", lw=0.5, ls="--", zorder=0)

    ax.set_xlabel("Fisher's Z", fontsize=10, labelpad=6)
    ax.set_title(title, fontsize=11, fontweight="bold", pad=10)

    # Legend annotation
    ax.annotate("● Control  → SCZ", xy=(0.98, 0.02),
                xycoords="axes fraction", ha="right", va="bottom",
                fontsize=9, color="#222222",
                bbox=dict(boxstyle="round,pad=0.3", fc="white",
                          ec="0.8", alpha=0.9))

    # Clean up spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)

    # Margin
    x_min = min(df["z_ctrl"].min(), df["z_scz"].min())
    x_max = max(df["z_ctrl"].max(), df["z_scz"].max())
    x_pad = (x_max - x_min) * 0.12 if x_max > x_min else 0.1
    ax.set_xlim(x_min - x_pad, x_max + x_pad)

    plt.tight_layout()
    if output_path:
        for ext in ['.pdf', '.png', '.svg']:
            plt.savefig(f"{output_path}{ext}", dpi=300, bbox_inches="tight")
    plt.show()
