import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from scipy import stats
import os
import re

def get_protein_coding_genes_from_gtf(gtf_path):
    """get protein coding gene ENSEMBL id and symbol"""
    gtf = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=["seqname","source","feature","start","end","score","strand","frame","attribute"],
        low_memory=False
    )

    gtf_genes = gtf[gtf["feature"] == "gene"].copy()

    def parse_attr(attr, key):
        for field in attr.split(";"):
            field = field.strip()
            if field.startswith(key):
                return field.split('"')[1]
        return None

    gtf_genes["gene_id"]   = gtf_genes["attribute"].apply(lambda x: parse_attr(x, "gene_id"))
    gtf_genes["gene_name"] = gtf_genes["attribute"].apply(lambda x: parse_attr(x, "gene_name"))
    gtf_genes["gene_type"] = gtf_genes["attribute"].apply(lambda x: parse_attr(x, "gene_biotype"))  # release 114 is gene_biotype

    protein_coding_ids   = set(gtf_genes.loc[gtf_genes["gene_type"]=="protein_coding","gene_id"])
    protein_coding_names = set(gtf_genes.loc[gtf_genes["gene_type"]=="protein_coding","gene_name"])

    print(f"Protein-coding IDs: {len(protein_coding_ids):,}, names: {len(protein_coding_names):,}")

    return protein_coding_ids, protein_coding_names

def get_genes_protein_coding(adata, protein_coding_names):
    """Filter for Protein Coding Genes"""
    keep_pc = adata.var_names.isin(protein_coding_names)
    adata = adata[:, keep_pc].copy()
    return adata

def get_genes_min_frac(adata, min_exp=0.0, min_frac=0.30):
    """Returns genes that are expressed in at least min_frac of samples"""
    X = adata.X
    n = X.shape[0]

    if sp.issparse(X):
        hits = (X > min_exp).getnnz(axis=0)
        frac = np.asarray(hits, dtype=float).ravel() / float(n)
    else:
        frac = (np.asarray(X) > min_exp).sum(axis=0).astype(float) / float(n)

    mask = frac >= float(min_frac)
    genes_min_frac = adata.var_names[mask].tolist()
    return genes_min_frac

def get_genes_highly_variable(adata, threshold=0.25):
    """Returns top threshold fraction of genes by variance"""
    def _col_var(X):
        if sp.issparse(X):
            X = X.tocsr()
            mean = np.asarray(X.mean(axis=0)).ravel()
            sq_mean = np.asarray(X.multiply(X).mean(axis=0)).ravel()
            var = sq_mean - mean**2
            return np.clip(var, 0, None)
        else:
            return np.var(np.asarray(X), axis=0)

    X = adata.X
    vars_ = _col_var(X)

    keep_top = float(threshold)
    keep_top = 0 if keep_top < 0 else keep_top
    q = np.quantile(vars_, keep_top)

    mask = vars_ >= q
    genes_highly_variable = adata.var_names[mask].tolist()
    return genes_highly_variable

def _to_2d_array(X):
    if sp.issparse(X):
        return X.toarray()
    arr = np.asarray(X)
    return arr if arr.ndim == 2 else np.atleast_2d(arr)

def _bh_fdr(pvals):
    """Benjamini-Hochberg FDR q-value"""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    q = p * n / ranks
    qo = q[order]
    qo = np.minimum.accumulate(qo[::-1])[::-1]
    q[order] = np.clip(qo, 0, 1)
    return q

def get_degs(ad1, ad2, deg_alpha=0.05, deg_lfc_thresh=0.0, deg_eps=1e-8, deg_direction="adata2_over_adata1"):
    genes_common = ad1.var_names.intersection(ad2.var_names)
    ad1 = ad1[:, ad1.var_names.isin(genes_common)].copy()
    ad2 = ad2[:, ad2.var_names.get_indexer(ad1.var_names)].copy()

    A1 = _to_2d_array(ad1.X)
    A2 = _to_2d_array(ad2.X)

    t_stat, pvals = stats.ttest_ind(A2, A1, axis=0, equal_var=False, nan_policy='propagate')
    pvals = np.nan_to_num(pvals, nan=1.0, posinf=1.0, neginf=1.0)

    mean1 = A1.mean(axis=0)
    mean2 = A2.mean(axis=0)

    if deg_direction == "adata2_over_adata1":
        log2fc = np.log2((mean2 + deg_eps) / (mean1 + deg_eps))
    elif deg_direction == "adata1_over_adata2":
        log2fc = np.log2((mean1 + deg_eps) / (mean2 + deg_eps))
    else:
        raise ValueError("deg_direction must be 'adata2_over_adata1' or 'adata1_over_adata2'")

    is_deg = (pvals <= deg_alpha) & (np.abs(log2fc) >= deg_lfc_thresh)

    return is_deg, genes_common

def get_filtered_pseudobulk(adata_ctrl, adata_scz, protein_coding_names):
    # 1. Protein Coding Gene Filter
    print("- Applying Protein Coding Gene Filter...")
    adata_ctrl = get_genes_protein_coding(adata_ctrl, protein_coding_names)
    adata_scz = get_genes_protein_coding(adata_scz, protein_coding_names)
    print(f"# Genes after Protein Coding Gene Filter:{adata_ctrl.n_vars}, {adata_scz.n_vars}")
    
    # 2. Minimum Expression Filter
    print("- Applying Minimum Expression Filter...")
    genes_min_frac_ctrl = get_genes_min_frac(adata_ctrl)
    genes_min_frac_scz = get_genes_min_frac(adata_scz)
    genes_min_frac_merged = set(genes_min_frac_ctrl) | set(genes_min_frac_scz)
    print(f"# Genes Min Expression CTRL:{len(genes_min_frac_ctrl)}, SCZ:{len(genes_min_frac_scz)}, Merged:{len(genes_min_frac_merged)}")
    adata_ctrl_min_filt = adata_ctrl[:, adata_ctrl.var_names.isin(genes_min_frac_merged)]
    adata_scz_min_filt = adata_scz[:, adata_scz.var_names.isin(genes_min_frac_merged)]
    print(f"# Genes after Minimum Expression Filter:{adata_ctrl_min_filt.n_vars}, {adata_scz_min_filt.n_vars}")

    # 3. HVG/DEG Filter
    print("- Applying HVG/DEG Filter...")
    genes_highly_variable_ctrl = get_genes_highly_variable(adata_ctrl_min_filt)
    genes_highly_variable_scz = get_genes_highly_variable(adata_scz_min_filt)
    genes_highly_variable_merged = set(genes_highly_variable_ctrl) | set(genes_highly_variable_scz)
    print(f"# Genes HVG/DEG CTRL_HVG:{len(genes_highly_variable_ctrl)}, SCZ_HVG:{len(genes_highly_variable_scz)}, Merged_HVG:{len(genes_highly_variable_merged)}")

    is_deg, genes_common = get_degs(adata_ctrl_min_filt,
        adata_scz_min_filt,
        deg_alpha=0.05,
        deg_lfc_thresh=0.0,
        deg_eps=1e-8,
        deg_direction="adata2_over_adata1",
    )
    genes_DEG = set(genes_common[is_deg])
    print(f"# Genes DEG:{len(genes_DEG)}")

    genes_hvg_deg_merged = genes_DEG | genes_highly_variable_merged
    print(f"-- # final filtered genes:{len(genes_hvg_deg_merged)} --")

    adata_ctrl_filtered = adata_ctrl[:, adata_ctrl.var_names.isin(genes_hvg_deg_merged)]
    adata_scz_filtered = adata_scz[:, adata_scz.var_names.isin(genes_hvg_deg_merged)]

    return adata_ctrl_filtered, adata_scz_filtered

def save_df_h5ad(adata, out_path):
    df = adata.to_df()
    df.to_csv(out_path + ".csv")
    adata.write_h5ad(out_path + ".h5ad")

def main(gtf_path, cell_types, input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    protein_coding_ids, protein_coding_names = get_protein_coding_genes_from_gtf(gtf_path)
    
    for ct in cell_types:
        print(f"\nProcessing Cell Type: {ct}")
        try:
            adata_ctrl = sc.read_h5ad(f'{input_dir}/pb_{ct}_control_avg.h5ad')    
            adata_scz = sc.read_h5ad(f'{input_dir}/pb_{ct}_scz_avg.h5ad')
        except FileNotFoundError as e:
            print(f"Skipping {ct}: {e}")
            continue

        adata_ctrl_filtered, adata_scz_filtered = get_filtered_pseudobulk(adata_ctrl, adata_scz, protein_coding_names)
        
        out_path_ctrl = f'{output_dir}/pb_{ct}_control_filtered'
        out_path_scz = f'{output_dir}/pb_{ct}_scz_filtered'

        save_df_h5ad(adata_ctrl_filtered, out_path_ctrl)
        save_df_h5ad(adata_scz_filtered, out_path_scz)
        print(f"Saved filtered data for {ct} to {output_dir}")

