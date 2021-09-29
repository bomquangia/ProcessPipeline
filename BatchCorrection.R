Preprocess <- function(data, params, logger, seed) {
  # CheckInteger <- function(x) {
  #   return(!grepl("[^[:digit:]]", format(x, trim = TRUE, scientific = FALSE)))
  # }

  FindVariableFeatures <- function(data) {
    if (params$n_variable_features > 0) {
      data <- Seurat::FindVariableFeatures(data, selection.method = "vst", nfeatures = params$n_variable_features, verbose = FALSE)
    }
    return(data)
  }

  CorrectHarmony <- function(data) {
    data <- FindVariableFeatures(data)
    data <- Seurat::ScaleData(data, features = GetVariableFeatures(data, params), verbose = FALSE)
    data <- RunPCA(data, params, logger, seed)
    data <- harmony::RunHarmony(data, "bioturing_batch", plot_convergence = FALSE, verbose = FALSE, epsilon.harmony = -Inf, max.iter.harmony = 30)
    return(data)
  }

  CorrectCca <- function(data) {
    log4r::info(logger, paste('Running CCA'))
    objects <- Seurat::SplitObject(data, split.by = "bioturing_batch")
    misc <- data@misc
    objects <- lapply(objects, function(obj) {
      obj <- Seurat::NormalizeData(obj, verbose = FALSE)
      obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst",
                                          nfeatures = params$n_variable_features, verbose = FALSE)
    })
    batches <- levels(data@meta.data$bioturing_batch)
    ncol.list <- sapply(objects, ncol)
    if (any(ncol.list < 30)) {
      stop("At least one batch < 30 cells")
    }
    options(future.globals.maxSize=Inf)  # Prevent limits of RAM per thread.
    anchors <- Seurat::FindIntegrationAnchors(
      object.list = objects,
      dims = 1:30,
      k.filter = min(200, ncol.list),
      verbose = (.Platform$OS.type == "unix")
    )
    data <- Seurat::IntegrateData(anchorset = anchors, dims = 1:30, verbose = FALSE)
    data@assays$RNA@meta.features$bioturing_name <- objects[[1]]@assays$RNA@meta.features$bioturing_name # Recover gene name
    data@misc <- misc
    data@meta.data$bioturing_batch <- factor(data@meta.data$bioturing_batch, levels=batches)
    return(data)
  }

  CorrectMnn <- function(data) {
    data <- FindVariableFeatures(data)
    mnn.arg <- lapply(unique(data@meta.data$bioturing_batch), function(i) {
      data@assays$RNA@data[GetVariableFeatures(data, params), data@meta.data$bioturing_batch == i]
    })
    mnn.arg$k <- 30
    mnn.arg$cos.norm <- FALSE
    corrected.data <- do.call(batchelor::fastMNN, mnn.arg)
    corrected.data <- SingleCellExperiment::reducedDims(corrected.data)@listData$corrected
    data@reductions[["mnn"]] <- Seurat:::CreateDimReducObject(
      embeddings = corrected.data, key="mnn")
    return(data)
  }

  CorrectBatchEffect <- function(data) {
    # logger$cat("Removing batch effect with", terms$correct[[arg$correct_method]])
    log4r::info(logger, paste("Removing batch effect with", params$correct_method))
    batch.invalid = table(data@meta.data$bioturing_batch) < 5
    if (any(batch.invalid)) {
      batch.name <- names(batch.invalid)[batch.invalid]
      stop("Too few cells (fewer than 5 cells) in batch(s) ", paste(batch.name, collapse=", "))
    }
    data <- switch(
      params$correct_method,
      `mnn` = CorrectMnn(data),
      `cca` = CorrectCca(data),
      `harmony` = CorrectHarmony(data)
    )
    return(data)
  }

  if (params$unit == "umi") {
    log4r::info(logger, paste("Log normalizing RNA data..."))
    data <- Seurat::NormalizeData(data, assay = "RNA")
  } else if (params$unit == "lognorm") {
    log4r::info(logger, paste("Skipped normalization"))
  }

  if ("ADT" %in% names(data) && params$unit == "umi") {
    log4r::info(logger, paste("Normalizing ADT data..."))
    data <- NormalizeADT(data)
  }

  if (params$correct_method == "none") {
    if (params$n_variable_features > 0) {
      log4r::info(logger, paste("Finding top genes"))
      data <- Seurat::FindVariableFeatures(data, selection.method = "vst",
        nfeatures = params$n_variable_features, verbose = FALSE)
    }
  } else {
    data <- suppressWarnings(CorrectBatchEffect(data))
  }
  return(data)
}

RunPCA <- function(data, params, logger, seed) {
  log4r::info(logger, paste("Running PCA..."))
  set.seed(seed)
  is.integrated <- FALSE
  if ("integrated" %in% names(data@assays)) {
    Seurat::DefaultAssay(data) <- "integrated"
    is.integrated <- TRUE
  }
  data@assays$RNA@data@x[is.na(data@assays$RNA@data@x)] <- 0  # Issue 1788 of Seurat
  data <- Seurat::ScaleData(data, features = GetVariableFeatures(data, params), verbose = FALSE)
  if (is.integrated) {
    data <- Seurat::RunPCA(data, features = rownames(data@assays$integrated), verbose = FALSE)
  } else {
    if (ncol(data) < 3) {
      log4r::info(logger, paste("Too few data points to run PCA. Skipped."))
      emb <- t(data@assays$RNA@scale.data)
      colnames(emb) <- paste0("pca_", seq(ncol(emb)))
      data@reductions$pca <- Seurat::CreateDimReducObject(emb, assay = "RNA", key = "pca_")
    } else {
      npcs <- max(min(ceiling(ncol(data) / 2), 50), 2)
      data <- Seurat::RunPCA(data, npcs=npcs, verbose=FALSE, features=GetVariableFeatures(data, params))
    }
  }
  return(data)
}

GetVariableFeatures <- function(data, params) {
  if ((params$n_variable_features > 0) && length(Seurat::VariableFeatures(data)) == 0) {
      data <- Seurat::FindVariableFeatures(data, selection.method = "vst", nfeatures = params$n_variable_features, verbose = FALSE)
  }
  return(if (params$n_variable_features == 0) rownames(data) else Seurat::VariableFeatures(data))
}

GetBatchCorrection <- function(study_id, file_logger = NULL) {
  
  if (is.null(file_logger)) {
    file_logger = log4r::logger()
  }
  log4r::info(file_logger, paste("Getting batch correction method for", study_id))
  info <- GetParams(study_id, arg)
  correct_method <- info$correct_method
  if (!is.null(correct_method)) {
    log4r::info(file_logger, paste("Use batch correction method given by params:", correct_method))
    return(correct_method)
  }

  log4r::info(file_logger, paste("Batch correction method for this study not given by params, checking number of batch..."))
  
  info <- RunDiagnostics(study_id, arg = arg)
  original_n_batch <- info$original_n_batch
  
  if (original_n_batch == 1) {
    log4r::info(file_logger, paste("This study only has 1 batch, setting batch correction to `none`"))
    return("none")
  }
  
  log4r::info(file_logger, paste("More than one batch, using `harmony` as batch correction method"))
  return("harmony")
}


