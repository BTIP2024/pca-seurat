# pca-seurat
This package performs principal component analysis (PCA) on the scaled Seurat object. 

## Installation
The package can be installed using
```
devtools::install_github("BTIP/pca-seurat")
```

## Example
The output of this function would be 2D and 3D plots.
```
pca_seurat("after_scaling.rds")
```

## scRNAseq processing workflow 

The standard scRNAseq processing workflow with the R package Seurat consists of seven (7) steps. The output of this package and function should be used as input for the scRNAseq processing pipeline. 

The following are the repositories of the packages for every step of the pipeline:
1. QC and filtering: [qualitycontrolseurat package](https://github.com/BTIP2024/quality-control-seurat)
2. Normalization: [qualitycontrolseurat package](https://github.com/BTIP2024/quality-control-seurat)
3. Identification of highly variable features: [selectionscalingseurat package](https://github.com/BTIP2024/selection-scaling-seurat)
4. Scaling: [selectionscalingseurat package](https://github.com/BTIP2024/selection-scaling-seurat)
5. Linear Dimensionality Reduction (PCA): [pcaseurat package](https://github.com/BTIP2024/pca-seurat)
6. Clustering: [nonlinearreduction package](https://github.com/BTIP2024/non-linear-reduction)
7. Non-linear dimensionality reduction (t-SNE and UMAP): [nonlinearreduction package](https://github.com/BTIP2024/non-linear-reduction)

An overview of the pipeline and its outputs can be observed below:
![](https://github.com/user-attachments/assets/f5520299-5f90-48e2-b0d0-7a5acf71ca08)
