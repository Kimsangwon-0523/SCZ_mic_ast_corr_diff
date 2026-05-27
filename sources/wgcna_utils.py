import pandas as pd
import gseapy as gp
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

def get_module_dict(ct, z_summary_path, module_csv_path, z_summary_thres=2):
    """
    Returns a dictionary of module names/numbers mapped to gene lists,
    filtered by Zsummary.pres threshold.
    """
    zsummary = pd.read_csv(z_summary_path, index_col=0)
    module_list = zsummary[zsummary['Zsummary.pres'] > z_summary_thres].index.tolist()
    module_list.sort()

    # Number mapping
    module_mapping = {m: i + 1 for i, m in enumerate(module_list)}

    # Load module file
    df_module = pd.read_csv(module_csv_path, index_col=0)
    color_dict = df_module.groupby("dynamicColors")["genes"].apply(list).to_dict()

    # color_dict -> number_dict conversion
    number_dict = {
        module_mapping[color]: genes
        for color, genes in color_dict.items()
        if color in module_mapping
    }

    # Background genes: all genes that were assigned a color (including grey/excluded ones)
    background_genes = set(df_module["genes"])

    for key in sorted(number_dict.keys()):
        print(f'Module {key}: {len(number_dict[key])} genes')

    return number_dict, background_genes, module_mapping

def run_go_enrichment(gene_list, ct, module_idx, output_dir, gene_sets=['GO_Biological_Process_2021']):
    """
    Runs Enrichr for a given gene list and saves results to CSV.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    ora = gp.enrichr(
        gene_list=gene_list,
        gene_sets=gene_sets,
        organism='Human',
        outdir=None,
        cutoff=0.05
    )

    ora_results_df = ora.results.sort_values('Adjusted P-value')
    csv_path = os.path.join(output_dir, f'{ct}_module_{module_idx}_GO_BP.csv')
    ora_results_df.to_csv(csv_path, index=False)
    
    return ora_results_df

def plot_go_results(df_ora_res, ct, module_idx, output_dir, n_top=10, color_map=plt.cm.Purples):
    """
    Refined GO enrichment bar plot (horizontal).
    """
    top = df_ora_res.nsmallest(n_top, "Adjusted P-value").copy()
    top["-log10(Padj)"] = -np.log10(top["Adjusted P-value"])

    # Parse overlap: "48/621" -> numerator
    if "Overlap" in top.columns:
        top["Overlap_count"] = top["Overlap"].str.split("/").str[0].astype(int)
        top["Overlap_total"] = top["Overlap"].str.split("/").str[1].astype(int)

    # Clean term names
    top["Term_clean"] = top["Term"].str.replace(r"\s*\(GO:\d+\)", "", regex=True)

    # Sort
    top = top.sort_values("-log10(Padj)", ascending=True)

    # Styling
    plt.rcParams.update({
        "font.size": 9,
        "axes.linewidth": 0.6,
    })

    fig, ax = plt.subplots(figsize=(6, 4))
    cmap = plt.cm.Reds

    norm = mcolors.Normalize(
        vmin=top["-log10(Padj)"].min() * 0.5,
        vmax=top["-log10(Padj)"].max(),
    )
    bar_colors = cmap(norm(top["-log10(Padj)"].values))

    bars = ax.barh(
        range(len(top)),
        top["-log10(Padj)"].values,
        color=bar_colors,
        edgecolor="black",
        linewidth=0.4,
        height=0.7,
        zorder=3,
    )

    # Overlap count annotation
    # if "Overlap_count" in top.columns:
    #     for i, (val, count, total) in enumerate(
    #         zip(top["-log10(Padj)"], top["Overlap_count"], top["Overlap_total"])
    #     ):
    #         ax.text(
    #             val + 0.15, i,
    #             f"{count}/{total}",
    #             va="center", ha="left",
    #             fontsize=7, color="#666666",
    #         )

    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["Term_clean"].values, fontsize=8)
    ax.set_xlabel(r"$-\log_{10}$(adjusted $p$-value)", fontsize=10)

    x_max = top["-log10(Padj)"].max()
    ax.set_xlim(0, x_max * 1.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.grid(True, linewidth=0.3, color="#E8E8E8", zorder=0)
    ax.set_axisbelow(True)

    plt.title(f"{ct.capitalize()} Module {module_idx} GO Enrichment", weight="bold")
    plt.tight_layout()
    
    out_img = os.path.join(output_dir, f'{ct}_module_{module_idx}_GO_BP.png')
    plt.savefig(out_img, dpi=300, bbox_inches="tight")
    plt.close()
    return out_img
