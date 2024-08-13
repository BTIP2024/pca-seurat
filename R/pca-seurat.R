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
   others <- for_pca
   # 3d plots
   library(plotly)
   
   new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
   names(new.cluster.ids) <- levels(for_pca)
   for_pca <- Seurat::RenameIdents(for_pca, new.cluster.ids)
   plot_pca <- Seurat::DimPlot(for_pca, reduction = "pca", label = FALSE) + ggplot2::xlab("PCA 1") + ggplot2::ylab("PCA 2") + ggplot2::theme(axis.title = element_text(size = 18)) + ggplot2::guides(colour = guide_legend(override.aes = list(size = 10)))
   
   ggplot_pca <- ggplotly(plot_pca)
   
   htmltools::save_html(ggplot_pca, file = "pca_labeled.html")
   
   image <- Seurat::VizDimLoadings(others, dims = 1:2, reduction = "pca")
   ggplot2::ggsave(image, file = "pca_results.png", width = 15, height = 10)
   
   image1 <- Seurat::DimPlot(others, reduction = "pca", label = TRUE)
   ggplot2::ggsave(image1, file = "pca_plot_unlabeled.png", width = 15, height = 10)
   
   image2 <- Seurat::DimHeatmap(others, dims = 1, cells = 500, balanced = TRUE)
   ggplot2::ggsave(image2, file = "heatmap.png", width = 12, height = 12)
   
   image3 <- Seurat::DimHeatmap(others, dims = 1:10, cells = 500, balanced = TRUE)
   ggplot2::ggsave(image3, file = "heatmap_multiple.png", width = 20, height = 20)
   
   image4 <- Seurat::ElbowPlot(others)
   ggplot2::ggsave(image4, file = "elbowplot.png", width = 10)
}