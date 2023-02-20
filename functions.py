import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt



path_mtx_folder='data/filtered_gene_bc_matrices/hg19/'
path_file_h5ad='data/filtered_gene_bc_matrices/pbmc3k.h5ad'

def read_mtx(path_mtx_folder):
    """" Function for reading 10x formatted (filtered and normalized) single cell data
    Input: path to the folder containing matrix.mtx, barcodes.tsv and genes.tsv files
    Output: AnnData object

    """
    # Add module for reading unstructured data
    adata = sc.read_10x_mtx(path_mtx_folder,            # the directory with the `.mtx` file
                            var_names='gene_symbols',   # use gene symbols for the variable names (variables-axis index)
                            cache=True)                 # write a cache file for faster subsequent reading
    print('Reading mtx data')
    return adata


def read_h5ad(path_file_h5ad):
    """" Function for reading analyzed AnnData object
    Input: path to the folder containing h5ad file
    Output: AnnData object containing results from previous analysis
    
    """
    adata=sc.read_h5ad(path_file_h5ad)
    print('Reading h5ad data')
    return adata

def analysis_clustering(adata, n_neighbors=10, n_pcs=40):
    """" Function for streamlined analysis and clustering of single cell data
    Input:  - AnnData object (with filtered and normalized counts)
            - n_neighbors=10: number of neighbors
            - n_pcs=40: number of PC to consider

    Output: AnnData object containing analysis results and ready for visualisation
    """
    ## Logarithmize the data (if needed, check if already done)
    try:                                
        adata.uns['log1p']
        print('The data is already logaritmized')
    except:
        print('No log1p')
        print('Logarithmizing the data')
        sc.pp.log1p(adata)
        print('Logarithmizing the data - Done')
    
    ## Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5) 
    print('Identifyed highly variable genes')

    # ## Set the .raw attribute of the AnnData
    # adata.raw = adata

    ## Filter for highly_variable genes
    adata = adata[:, adata.var.highly_variable]

    ## Regress out effects of total counts per cell mito genes. Scale the data to unit variance.
    try:    
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        print('Scaling data on pct_counts_mt')
    except Exception:
        print('Could not find mitochondrial counts')
        try:    
            sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mito'])
            print('Scaling data on pct_counts_mito')
        except:
            print('Could not find mitochondrial counts')

    ##PCA 
    sc.tl.pca(adata, svd_solver='arpack') 
    print('Principal componenet analysis done')     

    ## Compute the neighborhood graph
    sc.pp.neighbors(adata, n_neighbors, n_pcs)
    print(f"Neighborhood graph done with n_neighbors={n_neighbors}") 

    ## Compute UMAP
    sc.tl.umap(adata) 
    print('UMAP done')

    ## Leiden clustering
    sc.tl.leiden(adata)  
    print(f"Clustering done with {len(adata.obs['leiden'].unique())} clusters")

    ## Find differentially expressed genes in each group
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    print('Ranking genes for each cluster')

    ## Return annotated AnnData object  
    return adata

def plot_umap(adata, map, color=['leiden']):
    """" Function for plotting composite figure with clustering results
        Input:
            - AnnData Object with clustering results
            - color plan for plotting
        Output: composite figure
            - umap with legend on the side
            - umap with legend plotted on clusteres
    """
    plt.rcParams['figure.constrained_layout.use'] = True
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4), layout="constrained")
    #plt.subplots_adjust(right=0.2)
    ax1_dict = sc.pl.embedding(adata, map, color=color, ax=ax1, show=False)
    ax2_dict = sc.pl.embedding(adata, map, color=color, ax=ax2, frameon=False, show=False,
                                legend_loc='on data', legend_fontsize=10, legend_fontoutline=1)
    return fig

def plot_dotplot(adata):
    """" Function for plotting dotplot for differentially expressed genes
        Input: AnnData Object with clustering results
        Output: dotplot for 5 top differentially expressed genes
    """
    plt.rcParams['figure.constrained_layout.use'] = True
    fig, (ax2) = plt.subplots(1, 1, figsize=(11,4), layout="constrained")
    ax2_dict = sc.pl.rank_genes_groups_dotplot(adata, n_genes=4, ax=ax2, show=False)
    return fig

# def plot_umap(adata, color=['leiden']):
#     fig, ax1 = plt.subplots(1, 1, figsize=(4,4))
#     ax1_dict = sc.pl.embedding(adata, 'X_umap', color=color, ax=ax1, frameon=False, show=False,
#                                 legend_loc='on data', legend_fontsize=10, legend_fontoutline=1)
#     return fig

# def plot_umap_gene(adata, color=['NKG7']):
#     fig, (ax2) = plt.subplots(1, 1, figsize=(4,4), gridspec_kw={'wspace':0.9})
#     ax2_dict = sc.pl.embedding(adata, 'X_umap', color=color, cmap='Reds', ax=ax2, show=False)
#     return fig


def plot_umap_gene_violin(adata, map, color=['NKG7']):
    """" Function for plotting composite figure with expression of a particular gene in different clusteres
        Input:
            - AnnData Object with clustering results
            - gene name 
        Output: composite figure
            - umap plot: expression of a chosen gene
            - violin plot: expression of a chosen gene
    """
    plt.rcParams['figure.constrained_layout.use'] = True
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11,4), gridspec_kw={'width_ratios': [2,3]}, layout="constrained")
    ax1_dict = sc.pl.embedding(adata, map, color=color, cmap='Reds', ax=ax1, show=False)
    ax2_dict = sc.pl.violin(adata, color, groupby='leiden', frameon=False, ax=ax2, show=False)
    return fig

def diff_exp_table(adata):
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    table=pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']}).head(20)
    return table

def rename_clusters():

    pass

def sub_setting():

    pass

def save_data():

    pass


if __name__ == '__main__':
    print(f"Simple Cell Explorer")
