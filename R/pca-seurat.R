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
   for_pca <- Seurat::RunPCA(for_pca, features = Seurat::VariableFeatures(for_pca))
   others <- for_pca
   
   # 2D plots no labels
   library(plotly)
   
   new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
   names(new.cluster.ids) <- levels(for_pca)
   for_pca <- Seurat::RenameIdents(for_pca, new.cluster.ids)
   plot_pca <- Seurat::DimPlot(for_pca, reduction = "pca", label = FALSE) + ggplot2::xlab("PCA 1") + ggplot2::ylab("PCA 2") + ggplot2::theme(axis.title = element_text(size = 18)) + ggplot2::guides(colour = guide_legend(override.aes = list(size = 10)))
   
   ggplot_pca <- ggplotly(plot_pca)
   
   htmltools::save_html(ggplot_pca, file = "pca_unlabeled.html")
   
   #2D with labels
   for_pca <- Seurat::FindNeighbors(for_pca, dims=1:15)
   for_pca <- Seurat::FindClusters(for_pca, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
   with_labels1 <- Seurat::DimPlot(for_pca, group.by = "RNA_snn_res.0.3", label = TRUE)
   with_labels2 <- Seurat::DimPlot(for_pca, group.by = "RNA_snn_res.0.1", label = TRUE)
   
   htmltools::save_html(with_labels1, file = "pca03_labeled.html")
   htmltools::save_html(with_labels2, file = "pca01_labeled.html")
   
   # Other plots
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
