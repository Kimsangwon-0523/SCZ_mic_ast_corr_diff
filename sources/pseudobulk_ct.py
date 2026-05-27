# From donor level snRNA-seq data to donor level pseudobulk data

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import pandas as pd
import gc
from scipy.sparse import csr_matrix
import data_preprocess as dp
from tqdm import tqdm

def get_common_genes(path_list):
    # ---------- 1st PASS : get common genes ----------
    common_genes = None
    for p in tqdm(path_list, desc="Scanning genes"):
        adata = sc.read_h5ad(p, backed='r')      # read only metadata
        genes_here = set(adata.var_names)
        adata.file.close()                       # close the file
        if common_genes is None:
            common_genes = genes_here
        else:
            common_genes &= genes_here           # intersection

    if len(common_genes) == 0:
        raise ValueError("There are no common genes in all donors!")

    gene_names = sorted(common_genes)            # sorted order

    return gene_names

def get_pseudobulk_by_cell_type(gene_names, path_list, ct_col, cell_types):
    # ---------- 2nd PASS : get average expression of each cell type per donor ----------
    mean_bank = {ct: [] for ct in cell_types}
    donor_ids = []

    for p in tqdm(path_list, desc="Processing donors"):
        adata = sc.read_h5ad(p)                  # load all data
        adata = adata[:, gene_names]             # no KeyError
        donor_id = p.split("_")[-2:]
        donor_id = f"{donor_id[0][6:]}_{donor_id[1][:-5]}"
        donor_ids.append(donor_id)

        for ct in cell_types:
            mask = adata.obs[ct_col] == ct
            if mask.sum() == 0:
                raise ValueError(f"{p} has no '{ct}' cell.")
            mean_expr = adata[mask].X.mean(axis=0)
            if not isinstance(mean_expr, np.ndarray):
                mean_expr = np.asarray(mean_expr).ravel()
            mean_bank[ct].append(mean_expr)

        del adata
        gc.collect()

    # ---------- 3rd PASS : make AnnData object ----------
    # donor_ids = [f"donor_{i:03d}" for i in range(1, len(path_list)+1)]/
    pb_ct = {}
    for ct in cell_types:
        mat = np.stack(mean_bank[ct])            # (191, #genes_common)
        adata_ct = ad.AnnData(
            X   = mat,
            obs = pd.DataFrame(index=donor_ids),
            var = pd.DataFrame(index=gene_names)
        )
        pb_ct[ct] = adata_ct
    return pb_ct

def save_pb_ct(pb_ct, output_path, condition):
    """
    pb_ct: dict, {cell_type: AnnData}
    output_path: str, path to save
    """
    for ct in pb_ct.keys():
        adata = pb_ct[ct]
        adata.X = csr_matrix(adata.X)
        file_name = f"{output_path}/pb_{ct}_{condition}_avg.h5ad"
        adata.write_h5ad(file_name)
        print(f"Saved {ct} pseudobulk data to {file_name}")

def main(donor_metadata_path, donor_adata_path, pb_ct_path, ct_col, cell_types):
    # Ensure output directory exists
    os.makedirs(pb_ct_path, exist_ok=True)

    condition_dict = dp.parse_donor_metadata(donor_metadata_path)

    ctrl_path_list = [f"{donor_adata_path}/{sample_id}_unaffected.h5ad" for sample_id in condition_dict['Unaffected']]
    scz_path_list = [f"{donor_adata_path}/{sample_id}_affected.h5ad" for sample_id in condition_dict['Affected']]
    
    path_list = scz_path_list + ctrl_path_list
    gene_names = get_common_genes(path_list)

    pb_ct_ctrl = get_pseudobulk_by_cell_type(gene_names, ctrl_path_list, ct_col, cell_types)
    pb_ct_scz = get_pseudobulk_by_cell_type(gene_names, scz_path_list, ct_col, cell_types)

    save_pb_ct(pb_ct_ctrl, pb_ct_path, condition='control')
    print("Control pseudobulk data saved.")
    save_pb_ct(pb_ct_scz, pb_ct_path, condition='scz')
    print("Schizophrenia pseudobulk data saved.")

# if __name__ == "__main__":
#     donor_metadata_path = "/data/home/swkim0523/research/Schizophrenia_refine/data/raw/raw_data/SZvillage_donorMetadata.txt"
#     donor_adata_path = "/data/home/swkim0523/research/Schizophrenia_refine/data/snRNA_donor" 
#     pb_ct_path = "/data/home/swkim0523/research/Schizophrenia_refine/data/pb_ct" # output path
    
#     ct_col = 'cellclass'
#     cell_types = [
#         "gabaergic", "glutamatergic", "oligodendrocyte",
#         "microglia", "astrocyte", "polydendrocyte", "endothelia"
#     ]

#     main(donor_metadata_path, donor_adata_path, pb_ct_path, ct_col, cell_types)