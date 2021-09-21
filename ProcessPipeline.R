library(rhdf5)

NormalizeADT <- function(data) {
  data <- Seurat::NormalizeData(data, assay = "ADT", normalization.method = "CLR")
  data <- Seurat::ScaleData(data, assay = "ADT")
  return(data)
}

readMtxFromThaoH5 <- function(filepath) {
  ## TODO: Check file path exist

  counts <- Matrix::sparseMatrix(
    j = readThaoH5Slot(filepath, "/indices") + 1, # Why j, not i? How to differentiate from j and i?
    p = readThaoH5Slot(filepath, "/indptr"),
    x = readThaoH5Slot(filepath, "/data"),
    dims = readThaoH5Slot(filepath, "/shape"),
    dimnames = list(
      readThaoH5Slot(filepath, "/features"),
      readThaoH5Slot(filepath, "/barcodes")
    )
  )
  print(dim(counts))
  return(counts)
}

readThaoH5Slot <- function(filepath, slot) {
  return(rhdf5::h5read(filepath, slot, drop=TRUE))
}

GetMetadata <- function(filepath) {
  h5Info <-  rhdf5::h5ls(filepath)
  metaIdx <- h5Info$group == "/original_metadata"
  
  if (sum(metaIdx) == 0) {
    return(NULL)
  }
  metaName <- h5Info$name[metaIdx]
  print(metaName)
  meta_data <- lapply(paste0("/original_metadata/", metaName) , readThaoH5Slot, filepath=filepath) 
  names(meta_data) <- metaName
  return(meta_data)
}

GetPath <- function(study_id) {
  return(ConnectPath(arg$raw_path, paste0(study_id, '.hdf5')))
}

GetInfo <- function(study_id) {
  filepath <- GetPath(study_id)
  batch <- readThaoH5Slot(filepath, "/batch")
  species <- readThaoH5Slot(filepath, "/species")
  features <- readThaoH5Slot(filepath, '/features') 
  ADT_indices <- grepl("ADT-", features)
  return(list(batch=batch, species = species, ADT_indices=ADT_indices))
}

isCITESEQ <- function(study_id) {
  filepath <- GetPath(study_id)
  features <- readThaoH5Slot(filepath, '/features') 
  ADT_indices <- grepl("ADT-", features)
  hasADT <- sum(ADT_indices) > 0
  if (hasADT) {
    print(paste(study_id, "has ADT features"))
  }
  return(sum(ADT_indices) > 0)
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

GetBatchCorrection <- function(study_id) {
  info <- RunDiagnostics(study_id)
  correct_method <- info$correct_method
  original_n_batch <- info$original_n_batch
  if (original_n_batch == 1) {
    return("none")
  }
  if (is.null(correct_method) || correct_method %in% c("cca", "harmony", "mnn")) {
    return("harmony")
  }
  return("none")
}

SanityCheck <- function(filepath) {
  dims = readThaoH5Slot(filepath, "/shape")
  print("Dimensions")
  print(dims)
  
  features <- readThaoH5Slot(filepath, "/features")
  barcodes <- readThaoH5Slot(filepath, "/barcodes")
  print("Example features")
  n <- length(features)
  print(features[c(1:5, (n-5):n)])
  print("Example barcodes")
  print(barcodes[1:5])
  batch <- readThaoH5Slot(filepath, "/batch")
  print(paste0("Number of batches: ", length(unique(batch))))
  print(unique(batch))
}

GetStudyId <- function(filepath) {
  # filepath must be a path to a hdf5 file
  study_id <- sub('\\.hdf5$', '', basename(filepath)) 
  return(study_id)
}

CreateSeuratObj <- function(data, name) {
  feature_type <- rep("RNA", nrow(data))
  adt_idx <- grep("^ADT-", rownames(data))
  if (length(adt_idx) > 0) {
    adt_counts <- data[adt_idx, ]
    data <- data[-adt_idx, ] # remove adt expression
  }
  print(paste('Adding RNA assay. Dimension:', nrow(data), ncol(data)))

  obj <- Seurat::CreateSeuratObject(counts = data, project = name,
      min.cells = 0, min.features = 0)

  if (length(adt_idx) > 0) {
    print(paste('Adding ADT assay. Dimension:', nrow(adt_counts), ncol(adt_counts)))
    print("Removing ADT prefix...")
    rownames(adt_counts) <- gsub("ADT-", "", rownames(adt_counts))
    print("ADT feature names after removal")
    print(rownames(adt_counts)[1:5])
    obj[["ADT"]] <- Seurat::CreateAssayObject(counts = adt_counts[, colnames(obj)])
  }
  return(obj)
}

GetUnit <- function(study_id) {
  info <- RunDiagnostics(study_id)
  if (!is.null(info$unit)) {
    return(info$unit)
  }
  filepath <- GetPath(study_id)
  count_data <- readMtxFromThaoH5(filepath)
  browser()
  if (sum(abs(round(count_data@x[1:1000], 0) - count_data@x[1:1000])) > 0) { # Non integer
    return("not umi")
  }
  return("umi")  # Assuming that if count.data is integer, it is umi
}

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
    if (is.null(study.params[[key]])) {
      study.params[[key]] <- default.params[[key]]
    }
  }
  return(study.params)
}

RunPipeline <- function(study_id, arg) {

  filepath <- ConnectPath(arg$raw_path, paste0(study_id, '.hdf5'))
  output_path <- ConnectPath(arg$output.dir, paste0(study_id, '.bcs'))
  
  if (file.exists(output_path)) {
    print(paste("This study has already been process:", study_id))
    return(FALSE)
  }
  print(paste("Processing:", study_id))
  SanityCheck(filepath)

  params <- CombineParam(arg$default.params, arg$all.study.params[[study_id]])
  
  count_data <- readMtxFromThaoH5(filepath)
  
  if (any(duplicated(colnames(count_data)))) {
    print("There is duplicated barcodes")
    return(FALSE)
  }

  obj <- CreateSeuratObj(count_data, study_id)
  
  # TODO: Check this
  unit <- GetUnit(study_id)
  if ("ADT" %in% names(obj) && unit == "umi") {
    print("Normalizing ADT data...")
    obj <- NormalizeADT(obj)
  }

  # Handle Clonotype info

  meta_data <- GetMetadata(filepath)
  if (!is.null(meta_data)) {

    for (meta_name in names(meta_data)) {
      print(paste("Adding metadata:", meta_name))
      meta <- as.numeric(meta_data[[meta_name]])
      if (any(is.na(meta))) { # To avoid mistakenly forcing character into numeric
        meta <- meta_data[[meta_name]]
      }
      obj <- Seurat::AddMetaData(obj, meta, col.name = meta_name)
    }
  }

  obj <- Seurat::NormalizeData(obj)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures=min(params$n_variable_features, nrow(obj)))
  obj <- Seurat::ScaleData(obj)
  obj <- Seurat::RunPCA(obj)
  
  correct_method <- GetBatchCorrection(study_id)
  key <- 'pca'
  if (correct_method != "none" && params$correct_method %in% c("cca", "mnn", "harmony")) {
    print(paste("Run batch correction for:", study_id))
    batch <- readThaoH5Slot(filepath, "/batch")
    obj <- Seurat::AddMetaData(obj, batch, col.name = 'batch')
    obj <- harmony::RunHarmony(obj, "batch", plot_convergence = FALSE, verbose = TRUE, epsilon.harmony = -Inf, max.iter.harmony = 30)
    key <- 'harmony'
  } else {
    print(paste("Not correct batch for:", study_id))
  }
  
  mat <- obj@reductions[[key]]@cell.embeddings
  mat <- mat[, 1:min(ncol(mat), 50)]
  
  # Run TSNE and UMAP
  tsne.embeddings <- RunTSNE(mat, params)
  obj@reductions[['tSNE']] <- Seurat::CreateDimReducObject(embeddings=tsne.embeddings, assay="RNA", key=paste0(key, "_"))
  
  umap.embeddings <- RunUMAP(mat, params) 
  obj@reductions[['UMAP']] <- Seurat::CreateDimReducObject(embeddings=umap.embeddings, assay="RNA", key=paste0(key, "_"))
  
  # Run Graph-based cluster
  if (params$correct_method %in% c("cca", "mnn", "harmony") && params$correct_method %in% names(obj@reductions)) {
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
  rBCS::ExportSeurat(obj, output_path, overwrite=TRUE)
}

ProcessBatch <- function(study_list, output_dir, raw_path, old_data, arg) {
  study_list <- study_list
  output_dir <- output_dir
  raw_path <- raw_path
  old_data <- old_data
  lapply(study_list, RunPipeline, arg=arg)
}

GenerateImportArg <- function(study_id, bcs_dir = output_dir, write_dir = "/Users/bioturing/.BioTBData/Data/SingleCell/Study/") {
  
  info <- GetInfo(study_id)
  import_arg <- list()
  import_arg$group_name_aws <- NULL
  import_arg$input_path <- ConnectPath(bcs_dir, paste0(study_id, '.bcs'))
  import_arg$path_type <- "local"
  import_arg$batch_names <- NULL
  import_arg$input_id <- c(uuid::UUIDgenerate())
  import_arg$correct_method <- GetBatchCorrection(study_id)
  import_arg$type <- c("bcs")
  import_arg$unit <- "unknown" # TODO: Check this again 
  import_arg$email <- "vu@bioturing.com"
  import_arg$species <- info$species
  import_arg$log_path <- ConnectPath(write_dir, study_id, "tmp/submit.log")
  import_arg$dimred_method <- c("tsne", "umap")
  import_arg$dimred_perplexity <- 0
  import_arg$dimred_mode <- "slow"
  import_arg$quant <- list(method = "unknown", ref = "unknown")
  import_arg$seed <- 2409
  import_arg$platform <- "unknown"
  import_arg$filter <- list(cell=0, gene = c(0,0), mito = 100, top = 2000)
  import_arg$subcluster <- "none"
  import_arg$db_path <- "/Users/bioturing/.BioTBData/App/Model"
  import_arg$output_path <- ConnectPath(write_dir, study_id, "main")
  import_arg$create_adt_gallery <- TRUE
  import_arg$refIndex <- "unknown"
  import_arg$hash_id <- uuid::UUIDgenerate()
  import_arg$title <- study_id
  import_arg$tmpDir <-  ConnectPath(write_dir, study_id, "tmp")
  import_arg$zip_bin <- list(zip = "zip", unzip = "unzip")
  import_arg$version <- 16
  import_arg$bbrowser_version <- "2.10.23"
  import_arg$hash <- uuid::UUIDgenerate()
  return(import_arg)
}


# import_arg <- GenerateImportArg("GSE123022", bcs_dir = "/mnt2/vu/script/output/batch_4", write_dir='')

# lapply(study_list, function(study_id) {
#   import_arg <- GenerateImportArg(study_id, bcs_dir = "/mnt2/vu/script/output/batch_4", write_dir='')
#   jsonlite::write_json(import_arg, ConnectPath(output_dir, paste0(study_id, '.json')))
# })

# EXAMPLE CODE
# output_dir <- '/mnt2/vu/script/output'
# old_data <- "/mnt2/vu/script/old_data"
# raw_path <- "/mnt2/bioturing_data/raw_hdf5/"

# batch_4 = "/mnt2/bioturing_data/files/batch_4.tsv"
# study_list = as.character(read.table(file = batch_4, sep = '\t', header = FALSE)$V1)

# default.params <- jsonlite::fromJSON('/mnt2/vu/script/ProcessPipeline/default_params.json')
# all.study.params <- jsonlite::fromJSON('/mnt2/vu/script/ProcessPipeline/study_param.json')

# arg <- list(output.dir = output_dir, raw_path = raw_path, default.params = default.params, all.study.params = all.study.params)

# ProcessBatch(study_list, output_dir, raw_path, old_data, arg)