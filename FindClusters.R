#' Run graph-based and k-means clustering
#'
#' @param data a Seurat object
#' @param logger a logger object from a parent's function
FindClusters <- function(data, logger, arg) {
  log4r::info(logger, paste("Clustering..."))
  if (arg$correct_method %in% c("mnn", "harmony") && arg$correct_method %in% names(data@reductions)) {
    n.dim.use <- min(30, ncol(data@reductions[[arg$correct_method]]))
    data <- Seurat::FindNeighbors(data, dims = 1:n.dim.use, reduction = arg$correct_method)
  } else {
    n.dim.use <- min(30, ncol(data@reductions$pca))
    data <- Seurat::FindNeighbors(data, dims = 1:n.dim.use, reduction = "pca")
  }
  log4r::info(logger, paste("Running Louvain"))
  data <- Seurat::FindClusters(data, verbose = FALSE)
  data@meta.data$bioturing_graph <- as.numeric(data@meta.data$seurat_clusters)
  data@meta.data$bioturing_graph <- factor(
    paste("Cluster", data@meta.data$bioturing_graph),
    levels = paste("Cluster", sort(unique(data@meta.data$bioturing_graph)))
  )
  return(data)
}
