import anndata as ad
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import loompy as lp
from typing import Optional, List

def load_regulon_scores(loom_path: str) -> pd.DataFrame:
    """
    Load regulon scores from processed loom file.
    
    Args:
        loom_path (str): Path to processed loom file with AUCell scores
    Returns:
        pd.DataFrame: DataFrame containing regulon scores
    """
    with lp.connect(loom_path, mode="r", validate=False) as lf:
        auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
    return auc_mtx

def generate_embeddings(auc_mtx: pd.DataFrame, 
                       n_neighbors: int = 10,
                       metric: str = "correlation") -> ad.AnnData:
    """
    Generate UMAP and t-SNE embeddings from regulon scores.
    
    Args:
        auc_mtx (pd.DataFrame): Regulon scores matrix
        n_neighbors (int): Number of neighbors for embedding
        metric (str): Distance metric for neighbor computation
    Returns:
        ad.AnnData: Anndata object with computed embeddings
    """
    ad_auc_mtx = ad.AnnData(auc_mtx)
    sc.pp.neighbors(ad_auc_mtx, n_neighbors=n_neighbors, metric=metric)
    sc.tl.umap(ad_auc_mtx)
    sc.tl.tsne(ad_auc_mtx)
    return ad_auc_mtx

def plot_regulon_embeddings(adata: ad.AnnData,
                           embedding_key: str = "X_umap_aucell",
                           color_by: str = "Mutation",
                           save_path: Optional[str] = None) -> None:
    """
    Plot cell embeddings colored by metadata.
    
    Args:
        adata (ad.AnnData): Annotated data matrix
        embedding_key (str): Key for embedding coordinates
        color_by (str): Column name in adata.obs for coloring
        save_path (str, optional): Path to save the plot
    """
    sc.pl.embedding(adata, 
                   basis=embedding_key, 
                   color=color_by,
                   save=save_path is not None)
    if save_path:
        plt.savefig(save_path)
        plt.close()

def plot_regulon_heatmap(auc_mtx: pd.DataFrame,
                        groupby: str,
                        top_n_regulons: Optional[int] = None,
                        save_path: Optional[str] = None) -> None:
    """
    Generate clustered heatmap of regulon activities.
    
    Args:
        auc_mtx (pd.DataFrame): Regulon scores matrix
        groupby (str): Column name to group cells by
        top_n_regulons (int, optional): Number of top variable regulons to show
        save_path (str, optional): Path to save the plot
    """
    mean_auc = auc_mtx.groupby(groupby).mean()
    
    if top_n_regulons:
        var_scores = mean_auc.var(axis=0)
        top_regulons = var_scores.nlargest(top_n_regulons).index
        mean_auc = mean_auc[top_regulons]
    
    plt.figure(figsize=(15, 6.5))
    g = sns.clustermap(mean_auc,
                      cmap="Blues",
                      xticklabels=True,
                      yticklabels=True)
    
    if save_path:
        g.savefig(save_path)
        plt.close()

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Visualize pySCENIC results')
    parser.add_argument('--loom_path', required=True, help='Path to processed loom file')
    parser.add_argument('--adata_path', required=True, help='Path to original annotated data')
    parser.add_argument('--output_dir', required=True, help='Directory for output files')
    parser.add_argument('--group_by', default='Mutation', help='Column name for grouping cells')
    parser.add_argument('--top_n_regulons', type=int, default=50, help='Number of top regulons to show')
    
    args = parser.parse_args()
    
    # Load data
    auc_mtx = load_regulon_scores(args.loom_path)
    adata = sc.read(args.adata_path)
    
    # Generate embeddings
    ad_auc_mtx = generate_embeddings(auc_mtx)
    adata.obsm["X_umap_aucell"] = ad_auc_mtx.obsm["X_umap"]
    
    # Create visualizations
    os.makedirs(args.output_dir, exist_ok=True)
    
    # UMAP plot
    plot_regulon_embeddings(adata, 
                           color_by=args.group_by,
                           save_path=f"{args.output_dir}/umap_plot.pdf")
    
    # Regulon activity heatmap
    auc_mtx[args.group_by] = adata.obs[args.group_by]
    plot_regulon_heatmap(auc_mtx,
                        groupby=args.group_by,
                        top_n_regulons=args.top_n_regulons,
                        save_path=f"{args.output_dir}/regulon_heatmap.pdf")

if __name__ == "__main__":
    main()
