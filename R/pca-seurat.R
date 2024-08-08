#' PCA on Seurat Data
#' 
#' This package performs linear dimensional reduction on the seurat object, generating png files
#' 
#' @param input is the file after scaling
#' @examples
#' pca_seurat("after_scaling.rds")
#' @export

pca_seurat <- function(input){
   for_pca <- readRDS(input)
   for_pca <- Seurat::RunPCA(for_pca, features = Seurat::VariableFeatures(object = for_pca))
   print(for_pca[["pca"]], dims = 1:5, nfeatures = 5)
   
   image <- Seurat::VizDimLoadings(for_pca, dims = 1:2, reduction = "pca")
   ggplot2::ggsave(image, file = "pca.png", width = 15, height = 10)
   
   image2 <- Seurat::DimHeatmap(for_pca, dims = 1, cells = 500, balanced = TRUE)
   ggplot2::ggsave(image2, file = "heatmap.png", width = 12, height = 12)
   
   image3 <- Seurat::ElbowPlot(for_pca)
   ggplot2::ggsave(image3, file = "elbowplot.png", width = 10)
   
   image4 <- Seurat::FindNeighbors(for_pca, dims= 1:15)
   image4 <- Seurat::FindClusters(image4, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
   saveRDS(image4, file = "clustersperresolution.rds")
   
   image4 <- Seurat::DimPlot(image4, group.by = "RNA_snn_res.0.3", label = TRUE)
   
   ggplot2::ggsave(image4, file = "dimplot.png", width = 12, height = 10)
}
