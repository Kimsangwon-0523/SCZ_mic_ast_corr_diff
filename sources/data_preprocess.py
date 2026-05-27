## from batch snRNA-seq data, get adata of each sample.

import os
import h5py
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csc_matrix
import scipy.sparse as sp
from anndata import AnnData


def load_adata(file_path):
    """ From H5 format data, generate adata object."""
    with h5py.File(file_path, 'r') as f:
        barcodes = f['matrix/barcodes'][:]
        gene_names = f['matrix/features']['name'][:]  
        data = f['matrix/data'][:]
        indices = f['matrix/indices'][:]
        indptr = f['matrix/indptr'][:]
        shape = f['matrix/shape'][:]

    barcodes = [b.decode('utf-8') if isinstance(b, bytes) else b for b in barcodes]
    gene_names = [b.decode('utf-8') if isinstance(b, bytes) else b for b in gene_names]
    X = csc_matrix((data, indices, indptr), shape=shape)
    adata = sc.AnnData(X=X.T, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=gene_names))
    
    return adata

def merge_anndata(adatas):
    """ Function for merging adata objects."""
    merged = sc.concat(adatas, join='outer', label='rxn', keys=[f'rxn_{i}' for i in range(len(adatas))])
    return merged

def get_batch_adata(batch_name, raw_data_path):
    # batch_name: the name of the batch
    # raw_data_path: the path to the raw data
    # return adata object of each batch
    if batch_name == 'BA46_05-07-2019':
        filenames = [f'{raw_data_path}{batch_name}_rxn{i}_MyelinRemoved.h5' for i in range(2, 9)]
    else:
        filenames = [f'{raw_data_path}{batch_name}_rxn{i}.h5' for i in range(1, 9)]

    adata_list = []
    for file_name in filenames:
        adata_list.append(load_adata(file_name))

    merged = merge_anndata(adata_list)
    return merged

def filter_joint_metadata(joint_metadata, batch_name):
    # batch_name: the name of the batch
    # joint_metadata: the joint metadata
    # return filtered metadata
    joint_metadata_batch = pd.DataFrame()
    for i in range(1,9):
        joint_metadata_batch = pd.concat([joint_metadata_batch, joint_metadata[joint_metadata['PREFIX'] == f'{batch_name}_rxn{i}']])

    return joint_metadata_batch

def parse_donor_metadata(donor_metadata_path):
    """ Parse donor metadata"""
    metadata = pd.read_csv(donor_metadata_path, sep='\t')
    condition_dict = {'Affected':[], 'Unaffected':[]}
    for idx, row in metadata.iterrows():
        if row['Schizophrenia'] == 'Affected':
            condition_dict['Affected'].append(row['Donor'])
        elif row['Schizophrenia'] == 'Unaffected':
            condition_dict['Unaffected'].append(row['Donor'])
    
    return condition_dict

def save_adatas(donor_adata, donor_metadata_path, donor_adata_path):
    condition_dict = parse_donor_metadata(donor_metadata_path)

    for donor in donor_adata.keys():
        if donor in condition_dict['Unaffected']:
            donor_adata[donor].write_h5ad(f"{donor_adata_path}/{donor}_unaffected.h5ad")
        else:
            donor_adata[donor].write_h5ad(f"{donor_adata_path}/{donor}_affected.h5ad")


def sum_cells_by_barcode(adata):
    """If there are overlapping barcodes, sum them up."""
    # cell barcodes
    barcodes = adata.obs.index.values
    # unique barcode and inverse array
    unique_barcodes, inverse = np.unique(barcodes, return_inverse=True)
    
    if sp.issparse(adata.X):
        n_cells = len(barcodes)
        n_groups = len(unique_barcodes)
        rows = inverse          
        cols = np.arange(n_cells)   
        data = np.ones(n_cells)
        # (n_groups, n_cells)
        grouping_matrix = sp.coo_matrix((data, (rows, cols)), shape=(n_groups, n_cells))
        # multiply grouping matrix and adata.X, then sum them up by the group
        summed_X = grouping_matrix.dot(adata.X).tocsr()
    else:
        # if dense matrix, use pandas DataFrame's groupby to sum up
        df_X = pd.DataFrame(adata.X, index=adata.obs.index)
        summed_X_df = df_X.groupby(level=0).sum()
        summed_X = summed_X_df.values
        unique_barcodes = summed_X_df.index.values
    
    # the metadata for the same barcode group uses the first value
    new_obs = adata.obs.groupby(adata.obs.index).first().loc[unique_barcodes]
    
    return AnnData(X=summed_X, obs=new_obs, var=adata.var)


def main(donor_metadata_path, joint_metadata_path, donor_adata_path, raw_data_path):
    # Ensure the output directory exists
    os.makedirs(donor_adata_path, exist_ok=True)
    
    joint_metadata = pd.read_csv(joint_metadata_path, sep='\t')

    batches = {'batch1':'BA46_05-07-2019', 
           'batch2':'BA46_05-21-2021',
           'batch3':'BA46_05-28-2021',
           'batch4':'BA46_07-10-2019',
           'batch5':'BA46_08-09-2019',
           'batch6':'BA46_08-27-2019',
           'batch7':'BA46_09-02-2021',
           'batch8':'BA46_09-03-2021',
           'batch9':'BA46_09-18-2019',
           'batch10':'BA46_09-22-2021',
           'batch11':'BA46_10-16-2019',
           }
    
    for batch_num in batches.keys():
        adata_batch = get_batch_adata(batches[batch_num], raw_data_path)
        joint_metadata_batch = filter_joint_metadata(joint_metadata, batches[batch_num])
        cell_barcodes = joint_metadata_batch['CELL_BARCODE'].tolist()
        adata_filtered = adata_batch[adata_batch.obs.index.isin(cell_barcodes)].copy()
        adata_summed = sum_cells_by_barcode(adata_filtered)
        adata = adata_summed
        joint_metadata_unique = joint_metadata_batch.drop_duplicates(subset='CELL_BARCODE').set_index('CELL_BARCODE')
        adata.obs['DONOR'] = adata.obs.index.map(joint_metadata_unique['DONOR'])
        adata.obs['cellclass'] = adata.obs.index.map(joint_metadata_unique['cellclass'])
        adata.obs['subclass'] = adata.obs.index.map(joint_metadata_unique['subclass'])
        donor_categories = adata.obs['DONOR'].unique()
        donor_adata = {donor: adata[adata.obs['DONOR'] == donor].copy() for donor in donor_categories}
        
        save_adatas(donor_adata, donor_metadata_path, donor_adata_path)
        print(batch_num, " ended")

if __name__=="__main__":
    main()