import pandas as pd
import numpy as np
import re
from scipy.stats import fisher_exact, spearmanr, binomtest
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from tqdm import tqdm

def parse_gwas_genes(series):
    """Parse gene names from various delimiters in GWAS catalog."""
    genes = set()
    for val in series.dropna():
        parts = re.split(r'[,;\s]+|(?<!\w)-(?!\w)', str(val))
        for gene in parts:
            gene = gene.strip()
            if gene and gene not in ('NR', 'NA', '') and not gene.isdigit():
                genes.add(gene)
    return genes

def load_gwas_genes(gwas_path):
    """Load and parse schizophrenia GWAS genes."""
    df_gwas = pd.read_csv(gwas_path, sep='\t', low_memory=False)
    # Using reported genes as per legacy
    gwas_genes = parse_gwas_genes(df_gwas['REPORTED GENE(S)'])
    return gwas_genes

def calculate_gwas_enrichment(module_dict, gwas_genes, background_genes):
    """
    Perform Fisher's exact test for GWAS gene enrichment in each module.
    """
    # Legacy calculation: background is ONLY the genes in the analyzed modules
    all_module_genes = set().union(*module_dict.values())
    
    gwas_genes_filtered = set(gwas_genes).intersection(all_module_genes)
    total_background_count = len(all_module_genes)
    
    results = []
    for mod_id, genes in module_dict.items():
        mod_genes = set(genes).intersection(background_genes)
        a = len(mod_genes.intersection(gwas_genes_filtered)) # In Module & In GWAS
        b = len(mod_genes) - a                             # In Module & Not GWAS
        c = len(gwas_genes_filtered) - a                   # Not In Module & In GWAS
        d = total_background_count - a - b - c             # Not In Module & Not GWAS
        print(mod_id, a,b,c,d)
        
        oddsratio, pvalue = fisher_exact([[a, b], [c, d]], alternative='greater')
        
        results.append({
            'Module': mod_id,
            'Module_Gene_Count': len(mod_genes),
            'GWAS_Overlap_Count': a,
            'Overlap_Ratio': a / len(mod_genes) if len(mod_genes) > 0 else 0,
            'OddsRatio': oddsratio,
            'P_value': pvalue
        })
        
    df_fisher = pd.DataFrame(results)
    _, p_adj, _, _ = multipletests(df_fisher['P_value'], method='fdr_bh')
    df_fisher['P_adj'] = p_adj
    return df_fisher

def calculate_rho_series(gwas_or_series, delta_z_df):
    """
    Calculate Spearman rho for each target module (rows of delta_z_df)
    against GWAS ORs of reference modules (columns of delta_z_df).
    """
    rho_results = {}
    p_results = {}
    common_ref_mods = gwas_or_series.index.intersection(delta_z_df.columns)
    
    for target_mod in delta_z_df.index:
        dz_vals = delta_z_df.loc[target_mod, common_ref_mods]
        or_vals = gwas_or_series.loc[common_ref_mods]
        
        rho, pval = spearmanr(dz_vals.values, or_vals.values)
        rho_results[target_mod] = rho
        p_results[target_mod] = pval

    return pd.Series(rho_results), pd.Series(p_results)

def run_permutation_test(gwas_or_series, delta_z_df, n_perm=100, seed=42, alternative='greater'):
    """
    Shuffle GWAS ORs and calculate mean Spearman rho null distribution.
    """
    np.random.seed(seed)
    rho_obs_series, p_obs_series = calculate_rho_series(gwas_or_series, delta_z_df)
    if len(rho_obs_series) == 0:
        return 0, np.zeros(n_perm), 1.0
        
    t_obs = rho_obs_series.mean()
    
    t_null = []
    common_ref_mods = gwas_or_series.index.intersection(delta_z_df.columns)
    
    for _ in tqdm(range(n_perm), desc="Permutation"):
        shuffled_or = gwas_or_series.copy()
        shuffled_or.loc[common_ref_mods] = np.random.permutation(shuffled_or.loc[common_ref_mods].values)
        
        rhos_perm = []
        for target_mod in delta_z_df.index:
            dz_vals = delta_z_df.loc[target_mod, common_ref_mods]
            or_vals = shuffled_or.loc[common_ref_mods]
            
            mask = (~dz_vals.isna()) & (~or_vals.isna())
            if mask.sum() > 2:
                try:
                    rho, _ = spearmanr(or_vals[mask], dz_vals[mask])
                    if not np.isnan(rho):
                        rhos_perm.append(rho)
                except:
                    pass
        
        if len(rhos_perm) > 0:
            t_null.append(np.mean(rhos_perm))
        else:
            t_null.append(0)
        
    t_null = np.array(t_null)
    if alternative == 'greater':
        pval = np.mean(t_null >= t_obs)
    else:
        pval = np.mean(t_null <= t_obs)
    
    return t_obs, t_null, pval

def plot_gwas_rho_barplot(rho_series, t_obs, title, save_path=None):
    """Panel (b): Bar plot of Spearman rho per module."""
    plt.rcParams.update({
        'font.family': 'Liberation Serif',
        'font.size': 12,
        'axes.linewidth': 0.8,
    })
    
    # Color mapping
    cmap = cm.get_cmap('RdBu_r')
    vmax = max(abs(rho_series.min()), abs(rho_series.max())) if len(rho_series) > 0 else 1
    norm = mcolors.Normalize(vmin=-vmax, vmax=vmax)
    bar_colors = [cmap(norm(r)) for r in rho_series.values]
    
    fig, ax = plt.subplots(figsize=(5, 6))
    y_pos = np.arange(len(rho_series))
    ax.barh(y_pos, rho_series.values, color=bar_colors, edgecolor='white', height=0.7, linewidth=0.5)
    
    ax.axvline(0, color='black', linewidth=0.8)
    ax.axvline(t_obs, color='black', linewidth=1.5, linestyle='--', zorder=5)
    
    ax.annotate(
        f'Mean \u03c1 = {t_obs:.2f}',
        xy=(t_obs, len(rho_series)-0.5), xytext=(t_obs + 0.15, len(rho_series) + 0.5),
        fontsize=11, color='black', # fontweight='bold',
        arrowprops=dict(arrowstyle='->', color='black', lw=1.2, connectionstyle="arc3,rad=-0.2")
    )
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(rho_series.index, fontsize=12)
    ax.set_xlabel('Spearman \u03c1\n(Reference GWAS OR vs. \u0394Z)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.25, axis='x', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    if save_path:
        for ext in ['.pdf', '.png', '.svg']:
            plt.savefig(f"{save_path}{ext}", dpi=300, bbox_inches="tight")
    plt.show()

def plot_permutation_null(t_null, t_obs, pval, title, save_path=None):
    """Panel (c): Histogram of permutation null distribution."""
    plt.rcParams.update({
        'font.family': 'Liberation Serif',
        'font.size': 12,
    })
    
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(t_null, bins=60, color='#bdc3c7', edgecolor='white', linewidth=0.4, alpha=0.9, zorder=2)
    ax.axvline(t_obs, color='black', linewidth=2, linestyle='--', zorder=5)
    
    ax.set_xlabel('Null Distribution of Mean Spearman \u03c1', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.25, axis='y', linewidth=0.5)
    
    fig.canvas.draw()
    ylim = ax.get_ylim()
    ax.annotate(
        f'Mean \u03c1 = {t_obs:.3f}\nP-value = {pval:.3f}',
        xy=(t_obs, ylim[1] * 0.75),
        xytext=(t_obs + 0.02, ylim[1] * 0.85),
        fontsize=12, color='black', # #c0392b
        arrowprops=dict(arrowstyle='->', color='black', lw=1.2)
    )
    
    plt.tight_layout()
    if save_path:
        for ext in ['.pdf', '.png', '.svg']:
            plt.savefig(f"{save_path}{ext}", dpi=300, bbox_inches="tight")
    plt.show()
