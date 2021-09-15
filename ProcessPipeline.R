library(rhdf5)

readMtxFromThaoH5 <- function(filepath) {
  print(filepath)
  ## TODO: Check file path exist

  counts <- Matrix::sparseMatrix(
    j = rhdf5::h5read(filepath, "/indices", drop=TRUE) + 1, # Why j, not i? How to differentiate from j and i?
    p = rhdf5::h5read(filepath, "/indptr", drop=TRUE),
    x = rhdf5::h5read(filepath, "/data", drop=TRUE),
    dims = rhdf5::h5read(filepath, "/shape", drop=TRUE),
    dimnames = list(
      rhdf5::h5read(filepath, "/features", drop=TRUE),
      rhdf5::h5read(filepath, "/barcodes", drop=TRUE)
    )
  )  
  print(dim(counts))
  return(counts)
}

getInfo <- function(filepath) {
  batch <- rhdf5::h5read(filepath, "/batch", drop=TRUE)
  species <- rhdf5::h5read(filepath, "/species", drop=TRUE)
  return(list(batch=batch, species = species))
}

WriteSpMt <- function(filePath, groupName, mat) {
  if (!file.exists(filePath)) {
    rhdf5::h5createFile(filePath)
  }
  rhdf5::h5createGroup(filePath, groupName)
  rhdf5::h5createDataset(filePath, paste0(groupName, "/indices"), storage.mode = "integer",
                         dims = length(mat@i), chunk=min(10000, length(mat@i)), level=1)
  rhdf5::h5write(mat@i, filePath, paste0(groupName, "/indices"))
  rhdf5::h5createDataset(filePath, paste0(groupName, "/data"), storage.mode = "double",
                         dims = length(mat@x), chunk=min(10000, length(mat@x)), level=1)
  rhdf5::h5write(mat@x, filePath, paste0(groupName, "/data"))
  rhdf5::h5write(mat@p, filePath, paste0(groupName, "/indptr"))
  rhdf5::h5write(rownames(mat), filePath, paste0(groupName, "/features"))
  rhdf5::h5write(colnames(mat), filePath, paste0(groupName, "/barcodes"))
  rhdf5::h5write(dim(mat), filePath, paste0(groupName, "/shape"))
}

CreateStudyDir <- function(study.dir) {
  if (!dir.exists(study.dir)) {
    dir.create(study.dir)
  }
}

ConnectPath <- function(...) {
  res <- file.path(...)
  if (.Platform$OS.type == 'windows') {
    res <- gsub('\\/', '\\\\', res)
  }
  return(res)
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

SanityCheck <- function(filepath) {
  dims = rhdf5::h5read(filepath, "/shape", drop=TRUE)
  print("Dimensions")
  print(dims)
  
  features <- rhdf5::h5read(filepath, "/features", drop=TRUE)
  barcodes <- rhdf5::h5read(filepath, "/barcodes", drop=TRUE)
  print("Sample features")
  print(features[1:5])
  print("Sample barcodes")
  print(barcodes[1:5])
  batch <- rhdf5::h5read(filepath, "/batch", drop=TRUE)
  print(paste0("Number of batches: ", length(unique(batch))))
  print(unique(batch))
}

GetStudyId <- function(filepath) {
  # filepath must be a path to a hdf5 file
  study.id <- sub('\\.hdf5$', '', basename(filepath)) 
  return(study.id)
}

#CompareBarcodes <- function(obj, ) {
  #}
#AddOldMetadata <- function(old.study.path) {
  # compare barcodes.tsv
  
#}

CombineParam <- function(default.params, study.params) {
  # USAGE: Enable user to specify params to be used for a particular study, if differs from the global params
  # Example
  # global_param <- list(dims = 2, perplexity = 10, output.dir = output.dir, seed = 2409, correct_method='harmony')
  # study_param <- list(dims = 3, seed = 2409, correct_method='none')
  # CombineParam(global_param, study_param) --> list(dims = 3, seed = 2409, correct_method='none')

  if (is.null(study.params)) {
    study.params = list()
  }
  for (key in names(default.params)) {
    print(key)
    if (is.null(study.params[[key]])) {
      study.params[[key]] <- default.params[[key]]
    }
  }
  return(study.params)
}

RunPipeline <- function(study.id, arg) {
  filepath <- ConnectPath(arg$raw.path, paste0(study.id, '.hdf5'))
  output.path <- ConnectPath(arg$output.dir, paste0(study.id, '.bcs'))
  
  if (file.exists(output.path)) {
    print(paste("This study has already been process:", study.id))
    return(FALSE)
  }
  SanityCheck(filepath)
  # browser()
  params <- CombineParam(arg$default.params, arg$all.study.params[[study.id]])
  # Auto parse study id from filepath if not provided in param
  # study.id <- ifelse(!is.null(param$study.id), param.study.id, GetStudyId(filepath))
  # Check if `filepath`.bcs already exists
  
  count.data <- readMtxFromThaoH5(filepath)
  
  # If have batches --> remove batch effect with harmony
  study_info <- getInfo(filepath)
  hasMultiBatch <- unique(study_info$batch) > 1

  # Run PCA
  obj <- Seurat::CreateSeuratObject(count.data, min.cells=params$min.cells, min.features=params$min.features)

  obj <- Seurat::NormalizeData(obj)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures=min(params$n_variable_features, nrow(obj)))
  obj <- Seurat::ScaleData(obj)
  obj <- Seurat::RunPCA(obj)
  
  # If have batches --> remove batch effect with harmony
  study_info <- getInfo(filepath)
  hasMultiBatch <- length(unique(study_info$batch)) > 1
  key <- 'pca' # TODO: Handle case when key = 'harmony' for multibatch study
  if (hasMultiBatch && params$correct_method %in% c("mnn", "harmony")) { # Should check whether user want to run batch-correction for this particular study or not
    print(paste("This study has multibatch:", study.id))
    # browser()
    batch <- rhdf5::h5read(filepath, "/batch", drop=TRUE)
    obj <- Seurat::AddMetaData(obj, batch, col.name = 'batch')
    obj <- harmony::RunHarmony(obj, "batch", plot_convergence = FALSE, verbose = TRUE, epsilon.harmony = -Inf, max.iter.harmony = 30)
    key <- 'harmony'
  }
  
  mat <- obj@reductions[[key]]@cell.embeddings
  mat <- mat[, 1:min(ncol(mat), 50)]
  
  # Run TSNE and UMAP
  tsne.embeddings <- RunTSNE(mat, params)
  obj@reductions[['tSNE']] <- Seurat::CreateDimReducObject(embeddings=tsne.embeddings, assay="RNA", key=paste0(key, "_"))
  
  umap.embeddings <- RunUMAP(mat, params) 
  obj@reductions[['UMAP']] <- Seurat::CreateDimReducObject(embeddings=umap.embeddings, assay="RNA", key=paste0(key, "_"))
  
  # Run Graph-based cluster?
  if (params$correct_method %in% c("mnn", "harmony") && params$correct_method %in% names(obj@reductions)) {
    n.dim.use <- min(30, ncol(obj@reductions[[params$correct_method]]))
    obj <- Seurat::FindNeighbors(obj, dims = 1:n.dim.use, reduction = params$correct_method)
  } else {
    n.dim.use <- min(30, ncol(obj@reductions$pca))
    obj <- Seurat::FindNeighbors(obj, dims = 1:n.dim.use, reduction = "pca")
  }
  print("Running Louvain")
  obj <- Seurat::FindClusters(obj, verbose = FALSE)
  obj@meta.data$bioturing_graph <- as.numeric(obj@meta.data$seurat_clusters)
  obj@meta.data$bioturing_graph <- factor(
    paste("Cluster", obj@meta.data$bioturing_graph),
    levels = paste("Cluster", sort(unique(obj@meta.data$bioturing_graph)))
  )
  # TODO: Remove metadata: bioturing_graph and RNA_snn_res.0.8
  ## 
  # Import metadata
  # Check that new and old barcodes is the same, if yes, simply copy all metalist file
  # if not, this is gonna be mind-bending
  
  #
  rBCS::ExportSeuratObject(obj, output.path, overwrite=FALSE)
}

