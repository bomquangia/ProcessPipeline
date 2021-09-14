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
RunPipeline <- function(filepath, param) {
  count.data <- readMtxFromThaoH5(filepath)
  
  # Auto parse study id from filepath if not provided in param
  study.id <- ifelse(!is.null(param$study.id), param.study.id, GetStudyId(filepath))
  
  # If have batches --> remove batch effect with harmony
  study_info <- getInfo(filepath)
  hasMultiBatch <- unique(study_info$batch) > 1
  # Run PCA
  obj <- Seurat::CreateSeuratObject(count.data, min.cells=0, min.features=0)
  obj <- Seurat::NormalizeData(obj)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures=min(2000, nrow(obj)))
  obj <- Seurat::ScaleData(obj)
  obj <- Seurat::RunPCA(obj)
  obj <- Seurat::RunTSNE(obj, seed.use = 2409)
  
  
  key <- 'pca' # HTODO: Handle case when key = 'harmony' for multibatch study
  mat <- obj@reductions[[key]]@cell.embeddings
  mat <- mat[, 1:min(ncol(mat), 50)]
  
  # Run TSNE and UMAP
  tsne.embeddings <- RunTSNE(mat, param)
  obj@reductions[['tSNE']] <- Seurat::CreateDimReducObject(embeddings=tsne.embeddings, assay="RNA", key=paste0(key, "_"))
  
  umap.embeddings <- RunUMAP(mat, param) 
  obj@reductions[['UMAP']] <- Seurat::CreateDimReducObject(embeddings=umap.embeddings, assay="RNA", key=paste0(key, "_"))
  
  # Run Graph-based cluster?
  
  # Import metadata
  #Check that new and old barcodes is the same, if yes, simply copy all metalist file
  
  #if not, this is gonna be mind-bending
  #
  rBCS::ExportSeuratObject(obj, ConnectPath(param$output.dir, paste0(study.id, '.bcs')), overwrite=TRUE)
}