# Single_Cell_Explorer
Web application for analysis and exploration of single cell RNA-seq data

- Upload and read data (h5ad file)

- Analyze data
Basic Analysis workflow:
1. Find highly variable genes
2. Principal Component Analysis 
3. Compute the neighborhood graph (chose n_neighbors)
4. Compute UMAP
5. Clustering (method = Leiden) 
6. Find differentially expressed genes in each group (method = Wilcoxon)

- Visualize results
1. UMAP
2. UMAP expression of selected gene
3. Dotplot: 5 top differentially expressed genes for each cluster
4. Top 20 differentially expressed genes(csv table)

- Additional features: direct visualisation of analyzed h5ad files

![image](https://user-images.githubusercontent.com/116521950/217526448-86808783-2e62-47e0-a94a-caa1b4bb557b.png)




The app is simple but scalable - additional functionalities
- Additional analysis workflows
- Adaptation to read annotated files from Single Cell Atlas (EMBL-EBI)
- Subset analysis



