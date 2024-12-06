import os
import scanpy as sc
import numpy as np

def load_data(file_path):
    """
    Load single-cell RNA sequencing data from h5ad file.
    
    Args:
        file_path (str): Path to h5ad file
    Returns:
        AnnData: Loaded data object
    """
    return sc.read(file_path)

def preprocess_data(adata, cell_type=None, n_top_genes=2000):
    """
    Preprocess the data with normalization and filtering.
    
    Args:
        adata (AnnData): Input data
        cell_type (str): Cell type to subset, if any
        n_top_genes (int): Number of highly variable genes to select
    Returns:
        AnnData: Processed data
    """
    if cell_type:
        adata = adata[adata.obs['BroadCellType'] == cell_type, :].copy()
    
    # Normalize and transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Select highly variable genes
    sc.pp.highly_variable_genes(adata, 
                              flavor='seurat', 
                              n_top_genes=n_top_genes, 
                              subset=True)
    
    return adata

def filter_genes(adata, use_hvg=True, tf_list=None):
    """
    Filter genes based on variability and TF list.
    
    Args:
        adata (AnnData): Input data
        use_hvg (bool): Whether to use highly variable genes
        tf_list (list): List of transcription factors
    Returns:
        AnnData: Filtered data
    """
    if use_hvg:
        if tf_list is None:
            tf_list = []
        mask = (adata.var["highly_variable"] == True) | adata.var.index.isin(tf_list)
        adata = adata[:, mask]
    
    return adata

if __name__ == "__main__":
    # Configuration
    save_dir = '/path/to/your/data'
    input_file = os.path.join(save_dir, "aim1c_gex.h5ad")
    
    # Load and process data
    adata = load_data(input_file)
    adata = preprocess_data(adata, cell_type='DAN', n_top_genes=2000)
    adata = filter_genes(adata, use_hvg=True)
    
    # Save processed data
    adata.write_h5ad(os.path.join(save_dir, "processed_data.h5ad"))
