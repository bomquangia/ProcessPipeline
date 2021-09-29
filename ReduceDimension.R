ReduceDimension <- function(data, arg, logger, seed) {
  NeedPCA <- function(data) {
    if (arg$correct_method %in% c("harmony", "mnn")) {
      return(FALSE)
    } else {
      return(!"pca" %in% names(data@reductions))
    }
  }

  if (is.null(arg$dimred_mode)) {
    arg$dimred_mode <- "quick"
  }

  if (NeedPCA(data)) {
    data <- RunPCA(data, arg, logger, seed)
  }
  # if (!arg$dimred_method %in% names(data@reductions)) {
    # Small data
    if (ncol(data) <= 10) {
      # When data is too small to do anything
      if (!"pca" %in% names(data@reductions)) {
        log4r::info(logger, paste("Not enough data to perform dimensionality reduction"))
        return(data) # TODO: js should handle data with no dimred
      }

      # Use PC1 and PC2 from PCA as dimred
      log4r::info(logger, paste("Cannot run", arg$dimred_method, "with", ncol(data), "cell(s). Use PCA instead..."))
      emb <- data@reductions$pca@cell.embeddings[, 1:2]
      data@reductions$`PC1 and PC2` <- Seurat::CreateDimReducObject(emb)
      return(data)
    }

    log4r::info(logger, paste("Reducing dimensionality for RNA data with TSNE"))
    
    data <- ReduceDimensionFromObject(data, list(
      correct_method = arg$correct_method,
      method = "tsne",
      dims = 2,
      perplexity = arg$perplexity,
      seed = seed
    ))
    log4r::info(logger, paste("Reducing dimensionality for RNA data with UMAP"))
    data <- ReduceDimensionFromObject(data, list(
      correct_method = arg$correct_method,
      method = "umap",
      dims = 2,
      perplexity = arg$perplexity,
      seed = seed
    ))
  # }
  return(data)
}

ReduceDimensionFromObject <- function(obj, param) {
  key <- if (param$correct_method %in% c("harmony", "mnn")) param$correct_method else "pca"
  mat <- obj@reductions[[key]]@cell.embeddings
  mat <- mat[, 1:min(ncol(mat), 50)]
  mat <- if (param$method == "umap") RunUMAP(mat, param) else RunTSNE(mat, param)
  obj@reductions[[param$method]] <- Seurat::CreateDimReducObject(embeddings=mat, assay="RNA", key=paste0(key, "_"))
  return(obj)
}

RunTSNE <- function(data, param) {
  set.seed(param$seed)
  return(Rtsne::Rtsne(data, check_duplicates=FALSE, pca=FALSE,
                      dims=param$dims, perplexity=param$perplexity, verbose=FALSE)$Y)
}

RunUMAP <- function(data, param) {
  set.seed(param$seed)
  return(uwot::umap(data, n_neighbors = param$perplexity, n_components = param$dims))
}
