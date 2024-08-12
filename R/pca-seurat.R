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
   ggplot2::ggsave(image, file = "pca_results.png", width = 15, height = 10)
   
   image1 <- Seurat::DimPlot(for_pca, reduction = "pca") + ggplot2::theme(legend.position = "none")
   ggplot2::ggsave(image1, file = "pca_plot.png", width = 15, height = 10)
   
   image2 <- Seurat::DimHeatmap(for_pca, dims = 1, cells = 500, balanced = TRUE)
   ggplot2::ggsave(image2, file = "heatmap.png", width = 12, height = 12)
   
   image3 <- Seurat::ElbowPlot(for_pca)
   ggplot2::ggsave(image3, file = "elbowplot.png", width = 10)
   
}
