library(rhdf5)

#' Write annotations
#'
#' @param arg argument from nodeJS
#' @export
ImportMetadata <- function(arg.json) {
  MapBarcode <- function(barcodes, meta) {
    # Filter cells that are not in the dataset
    index.matches <- NoraSC::MatchBarcodes(meta[[1]], barcodes)
    
    meta <- meta[which(!is.na(index.matches)), , drop = FALSE]
    index.matches <- index.matches[!is.na(index.matches)]
    rownames(meta) <- NULL
    if (length(index.matches) == 0) {
      stop("Could not find any barcode that can match your data")
    }

    # Extend metadata table to fit the actual cell number
    is.numeric <- sapply(meta, IsNumeric)
    meta2 <- lapply(seq_along(is.numeric), function(i) {
      if (is.numeric[i]) {
        return(rep(NA, length(barcodes)))
      } else {
        return(rep("Unassigned", length(barcodes)))
      }
    })
    names(meta2) <- colnames(meta)
    meta2$barcodes <- barcodes
    for (i in seq_along(is.numeric)[-1]) {
      if (is.numeric[i]) {
        meta2[[i]][index.matches] <- as.numeric(meta[[i]])
      } else {
        meta2[[i]][index.matches] <- as.character(meta[[i]])
      }
    }
    meta2 <- as.data.frame(meta2, check.names=FALSE, stringsAsFactors=FALSE)

    return(meta2)
  }

  PrettifyNames <- function(x) {
    x[1] <- "barcodes"
    return(gsub("[\n|\r]", " ", x))
  }

  ReadMetadataTable <- function(filepath) {
    special.names <- c("Graph-based clusters", "Batch", "Total count", "Total expressed feature")
    data <- NoraSC::ReadPlainTextTable(filepath)
    colnames(data) <- PrettifyNames(colnames(data))
    types <- as.character(sapply(data, function(x) class(x)[1]))
    
    if (any(duplicated(data[[1]]))) {
      stop("Duplicated barcodes in ", basename(filepath))
    }
    valid.column <- sapply(seq_along(types), function(i) {
      if (types[i] == "numeric") {
        return(TRUE)
      } else {
        return(length(unique(data[[i]])) <= arg$unique_limit)
      }
    })
    valid.column[1] <- TRUE # always needs barcodes
    data <- data[, valid.column, drop=FALSE]
    data <- data[, !(names(data) %in% special.names), drop=FALSE]
    return(data)
  }

  CreateMetaMultiBatch <- function(input.dir) {
    meta.tmp <- lapply(1:length(input.dir), function(index) {
      meta.file <- input.dir[index]
      if (!is.na(meta.file)) {
        data <- ReadMetadataTable(meta.file)
        
        data[[1]] <- paste(index, data[[1]], sep = "_")
        return(data)
      }
      return(NULL)
    })
    meta.tmp <- meta.tmp[!sapply(meta.tmp, is.null)]
    meta <- plyr::rbind.fill(meta.tmp)
    return(meta)
  }

  CreateMeta <- function(input.dir, status) {
    if (status == "single") {
      
      meta <- ReadMetadataTable(input.dir[1])
    } else if (length(input.dir) == 1) {
      
      meta <- ReadMetadataTable(input.dir)
    } else {
      
      meta <- CreateMetaMultiBatch(input.dir)
    }
    if (is.null(ncol(meta)) || ncol(meta) < 2) {
      stop("No valid metadata is found")
    }
    return(meta)
  }

  DoProcess <- function() {
    barcodes <- as.character(rhdf5::h5read(ConnectPath(arg$data_path, "matrix.hdf5"), "bioturing/barcodes"))
    meta <- CreateMeta(arg$input_path, arg$type)
    
    omit_columns <- grep("Barcodes", colnames(meta), ignore.case = TRUE)
    meta <- cbind(meta[1], meta[ , -omit_columns])

    meta <- MapBarcode(barcodes, meta)
    # Create metalist (no IDs yet)
    meta <- lapply(2:ncol(meta), function(i) CreateMetadataObject(meta[[i]], colnames(meta)[i]))

    # Replace old metadata
    meta.old <- ReadJSON(ConnectPath(arg$data_path, "metadata", "metalist.json"))$content
    names <- sapply(meta, function(x) x$name)
    names.old <- sapply(meta.old, function(x) x$name)
    report <- list(replaced=names[names %in% names.old], new=names[!names %in% names.old])
    for (info in meta[names %in% names.old]) {
      i <- match(info$name, names.old)
      info.old <- meta.old[[i]]
      if (is.null(names(info.old$history))) {
        info$history <- c(info.old$history, CreateHistory(list(email=arg$email, description="Replaced from file")))
      } else {
        info$history <- CreateHistory(list(email=arg$email, description="Replaced from file"))
      }
      meta.old[[i]] <- info
    }

    # Add new metadata
    meta <- meta[!names %in% names.old]
    names(meta) <- sapply(seq_along(meta), function(x) GenerateUUID())
    meta <- lapply(meta, function(info) {
      info$history <- CreateHistory(list(email=arg$email, description="Import from file"))
      return(info)
    })
    meta.old <- c(meta.old, meta)

    # Write metalist
    CreateDir(arg$output_path)
    ids <- names(meta.old)
    meta.old <- lapply(ids, function(id) {
      info <- meta.old[[id]]

      # Force list before writing json
      info$type <- if (is.null(info$type)) "category" else info$type
      if (info$type == "category") {
        info$clusterName <- as.list(info$clusterName)
        info$clusterLength <- as.list(info$clusterLength)
      }

      if (info$name %in% report$replaced || info$name %in% report$new) {
        WriteJSON(info, ConnectPath(arg$output_path, paste0(id, ".json")), encrypt=IsBioTuring(arg$email), auto_unbox=TRUE)
      }
      info$clusters <- NULL
      return(info)
    })
    names(meta.old) <- ids
    WriteJSON(list(content=meta.old, version=1),
        ConnectPath(arg$output_path, "metalist.json"), auto_unbox=TRUE)

    return(jsonlite::toJSON(report))
  }
  arg <- ReadJSON(arg.json)
  
  return(try(DoProcess()))
}

CreateMetadataObject <- function(x, name) {
  info <- list(clusters=x, name=name)
  if (is(x, "numeric")) {
    info$type <- "numeric"
    return(info)
  } else {
    info$type <- "category"
    info$clusters[is.na(info$clusters)] <- "Unassigned" 
    info$clusterName <- names(sort(table(info$clusters), decreasing=TRUE)) # Sort groups by size
    info$clusterName <- unique(c("Unassigned", info$clusterName)) # Unassigned is always 1st
    info$clusters <- factor(info$clusters, levels=info$clusterName) # Factorize labels to indices
    info$clusters <- as.numeric(info$clusters) - 1 # Convert to JS base-0 indices
    info$clusterLength <- sapply(seq_along(info$clusterName), function(i) {
      return(sum(info$clusters == (i - 1)))
    })
  }
  return(info)
}


NormalizeADT <- function(data) {
  data <- Seurat::NormalizeData(data, assay = "ADT", normalization.method = "CLR")
  data <- Seurat::ScaleData(data, assay = "ADT")
  return(data)
}

IsEnsemblID <- function(x) {
  return(any(grepl("^ENS", x[1: min(50, length(x))])))
}

ReadMtxFromBioTuringRawH5 <- function(filepath) {
  ## TODO: Check file path exist

  counts <- Matrix::sparseMatrix(
    j = ReadBioTuringRawH5Slot(filepath, "/indices") + 1, # Why j, not i? How to differentiate from j and i?
    p = ReadBioTuringRawH5Slot(filepath, "/indptr"),
    x = ReadBioTuringRawH5Slot(filepath, "/data"),
    dims = ReadBioTuringRawH5Slot(filepath, "/shape"),
    dimnames = list(
      ReadBioTuringRawH5Slot(filepath, "/features"),
      ReadBioTuringRawH5Slot(filepath, "/barcodes")
    )
  )
  print(dim(counts))
  return(counts)
}

ReadBioTuringRawH5Slot <- function(filepath, slot) {
  return(rhdf5::h5read(filepath, slot, drop=TRUE))
}

GetMetadata <- function(study_id) {
  filepath <- GetPath(study_id)
  h5Info <-  rhdf5::h5ls(filepath)
  metaIdx <- h5Info$group == "/original_metadata"
  
  if (sum(metaIdx) == 0) {
    return(NULL)
  }
  metaName <- h5Info$name[metaIdx]
  meta_data <- lapply(paste0("/original_metadata/", metaName) , ReadBioTuringRawH5Slot, filepath=filepath) 
  names(meta_data) <- metaName
  
  return(meta_data)
}

AddMetadataFromFile <- function(arg, output_path, study_id) {
  unzip(output_path, exdir = arg$output_dir)
  setwd(arg$output_dir)
  old_dir <- unzip(output_path, list=TRUE)$Name[1]
  import_arg <- list(input_path=ConnectPath(arg$meta_dir, paste0(study_id, '.tsv')),
                  data_path=ConnectPath(old_dir, "/main"),
                  # raw_h5 = GetPath(study_id),
                  output_path= ConnectPath(arg$output_dir, old_dir, 'main/metadata'),
                  type = "single",
                  email = "vu@bioturing.com",
                  unique_limit = 100,
                  bbrowser_version = "2.9.23")
  ImportMetadata(jsonlite::toJSON(import_arg))
  zip::zip(paste0(study_id, '.bcs'), old_dir, compression_level=1)
  unlink(old_dir, recursive = TRUE)
}


GetPath <- function(study_id) {
  return(ConnectPath(arg$raw_path, paste0(study_id, '.hdf5')))
}

GetInfo <- function(study_id) {
  filepath <- GetPath(study_id)
  batch <- ReadBioTuringRawH5Slot(filepath, "/batch")
  species <- ReadBioTuringRawH5Slot(filepath, "/species")
  features <- ReadBioTuringRawH5Slot(filepath, '/features') 
  ADT_indices <- grepl("ADT-", features)
  return(list(batch=batch, species = species, ADT_indices=ADT_indices))
}

isCITESEQ <- function(study_id) {
  filepath <- GetPath(study_id)
  features <- ReadBioTuringRawH5Slot(filepath, '/features') 
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

GetBatchCorrection <- function(study_id, file_logger = NULL) {
  if (is.null(file_logger)) {
    file_logger = log4r::logger()
  }
  log4r::info(file_logger, paste("Getting batch correction method for", study_id))
  info <- CombineParam(arg$default_params, arg$all_study_params[[study_id]])
  # info <- arg$params[[study_id]]
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

SanityCheck <- function(filepath, file_logger=NULL) {
  if (is.null(file_logger)) {
    file_logger = log4r::logger()
  }
  dims = ReadBioTuringRawH5Slot(filepath, "/shape")
  # print("Dimensions")
  # print(dims)

  # log4r::info(file_logger, paste("Dimensions", dims[1], dims[2]))  

  features <- ReadBioTuringRawH5Slot(filepath, "/features")
  barcodes <- ReadBioTuringRawH5Slot(filepath, "/barcodes")

  n <- length(features)
  print(features[c(1:5, (n-5):n)])
  
  # print("Example features")
  log4r::info(file_logger, paste("Example features: "))
  log4r::info(file_logger, paste(do.call(paste, as.list(features[c(1:5, (n-5):n)]))))
  
  # print("Example barcodes")
  # print(barcodes[1:5])
  log4r::info(file_logger, paste("Example barcodes: "))
  log4r::info(file_logger, paste(do.call(paste, as.list(barcodes[1:5]))))

  batch <- ReadBioTuringRawH5Slot(filepath, "/batch")
  # print(paste0("Number of batches: ", length(unique(batch))))
  # print(unique(batch))

  log4r::info(file_logger, paste("Number of batches: ", length(unique(batch))))
  log4r::info(file_logger, paste("Batch names: "))
  log4r::info(file_logger, paste(do.call(paste, as.list(unique(batch)))))

}

CreateSeuratObj <- function(data, name, file_logger=NULL) {
  if (is.null(file_logger)) {
    file_logger = log4r::logger()
  }
  feature_type <- rep("RNA", nrow(data))
  adt_idx <- grep("^ADT-", rownames(data))
  if (length(adt_idx) > 0) {
    adt_counts <- data[adt_idx, ]
    data <- data[-adt_idx, ] # remove adt expression
  }
  # print(paste('Adding RNA assay. Dimension:', nrow(data), ncol(data)))
  log4r::info(file_logger, paste('Adding RNA assay. Dimension:', nrow(data), ncol(data)))
  obj <- Seurat::CreateSeuratObject(counts = data, project = name,
      min.cells = 0, min.features = 0)

  if (length(adt_idx) > 0) {
    
    log4r::info(file_logger, paste('Adding ADT assay. Dimension:', nrow(adt_counts), ncol(adt_counts)))
    log4r::warn(file_logger, paste("Removing `ADT-` prefix to avoid conflicting with rBCS..."))
    log4r::info(file_logger, paste("ADT feature names before removal"))
    log4r::info(file_logger, paste(do.call(paste, as.list(rownames(adt_counts)[1:5]))))
    rownames(adt_counts) <- gsub("^ADT-", "", rownames(adt_counts))
    log4r::info(file_logger, paste("ADT feature names after removal"))
    log4r::info(file_logger, paste(do.call(paste, as.list(rownames(adt_counts)[1:5]))))

    obj[["ADT"]] <- Seurat::CreateAssayObject(counts = adt_counts[, colnames(obj)])
  }
  return(obj)
}

GetClonotypeData <- function(study_id) {
  filepath <- GetPath(study_id)
  if (!'clonotype' %in% rhdf5::h5ls(filepath)$name) {
    return(NULL)
  }
  v_gene <- ReadBioTuringRawH5Slot(filepath, "/clonotype/v_gene")
  j_gene <- ReadBioTuringRawH5Slot(filepath, "/clonotype/j_gene")
  cdr3 <- ReadBioTuringRawH5Slot(filepath, "/clonotype/cdr3")
  chain <- ReadBioTuringRawH5Slot(filepath, "/clonotype/chain")
  barcode <- ReadBioTuringRawH5Slot(filepath, "/clonotype/barcode")
  raw_clonotype_id <- ReadBioTuringRawH5Slot(filepath, "/clonotype/raw_clonotype_id")
  full_length <- ReadBioTuringRawH5Slot(filepath, "/clonotype/full_length")
  productive <- ReadBioTuringRawH5Slot(filepath, "/clonotype/productive")
  # What about 'c_gene, d_gene, cdr3_nt, contig_id, is_cell, raw_concesus_id?
  return(as.data.frame(list(v_gene=v_gene, j_gene=j_gene, cdr3=cdr3, chain=chain, barcode=barcode, raw_clonotype_id=raw_clonotype_id, full_length=full_length, productive=productive)))
  # return(list(v_gene=v_gene, j_gene=j_gene, cdr3=cdr3, chain=chain, barcode=barcode, raw_clonotype_id=raw_clonotype_id, full_length=full_length, productive=productive))
}

FilterClonotype <- function(data, bc) {
  
  ### Filter by full_length and productive
  data <- data[
      as.logical((tolower(data$full_length) == 'true') * (tolower(data$productive) == 'true')), ]
  ### Filter by barcodes
  data$index <- match(data$barcode, bc) - 1  # User base 0 for JS
  data <- data[!is.na(data$index), ]
  ### Filter by cdr3
  data <- data[tolower(data$cdr3) != 'none', ]
  ### Validate
  if (nrow(data) == 0) {
    stop('No valid clonotype is found')
  }
  data$barcode <- NULL
  return(data)
}

CountFromTCRData <- function(data) {
  
  QueryAntigenInformation <- function(cdr3) {
    ag.info <- vdj[cdr3]
    mhc <- c()
    antigen <- c()
    for (i in seq(length(ag.info))) {
      info <- ag.info[[i]]
      if (is.null(info)) {
        mhc[[i]] <- 'None'
        antigen[[i]] <- 'None'
      } else {
        mhc[[i]] <- paste0(info$`MHC class`[1], ': ', info$`MHC A`[1],
            ' (A) - ', info$`MHC B`[1], ' (B) (', info$Species[1], ')')
        antigen[[i]] <- paste0(info$Epitope[1], ' (gene ',
            info$`Epitope gene`[1], ') (', info$`Epitope species`[1], ')')
      }
    }
    return(list(mhc = mhc, antigen = antigen))
  }

  GetChainInfo <- function(id.list, key) {
    return(sapply(id.list,
        function(x) paste(chain.count[[key]][as.numeric(x)], collapse = '<br>')))
  }

  GetPairName <- function(id.list) {
    return(sapply(id.list, function(x) {
      return(paste(
        sort(paste0(chain.count$chain[as.numeric(x)], ': ', chain.count$cdr3[as.numeric(x)])),
        collapse = '<br>'
      ))
    }))
  }

  ### Count by chain types
  print('Counting by chains')
  
  chain.count <- aggregate(data[, 'index'],
    data[, c('v_gene', 'j_gene', 'chain', 'cdr3')], sort)
  chain.count$count <- sapply(chain.count$x, length)
  chain.count <- chain.count[order(chain.count$count, decreasing = TRUE), ]
  chain.count$id <- do.call(paste, chain.count[, c('v_gene', 'j_gene', 'chain', 'cdr3')])
  chain.count$name <- paste0(chain.count$chain, ': ', chain.count$cdr3)

  ### Query VDJ database
  # print('Searching for antigen information')
  # ag.info <- QueryAntigenInformation(chain.count$cdr3)
  # chain.count$antigen <- ag.info$antigen
  # chain.count$mhc <- ag.info$mhc

  ### Index map
  print('Indexing chains')
  chain.map <- as.list(seq(nrow(chain.count)))
  names(chain.map) <- chain.count$id
  data$id <- unlist(chain.map[do.call(paste, data[, c('v_gene', 'j_gene', 'chain', 'cdr3')])])

  ### Count by pairs
  print('Counting by pair')
  
  clono.count <- aggregate(data[, 'id'], as.data.frame(data[, 'index']),
      function(x) paste(sort(x), collapse = ';'))
  clono.count <- aggregate(clono.count[, 'index', drop = FALSE],
      clono.count[, 'id', drop = FALSE], sort)
  clono.count$count <- sapply(clono.count$index, length)
  clono.count <- clono.count[order(clono.count$count, decreasing = TRUE), ]
  clono.count$id <- strsplit(clono.count$id, ';')

  ### Map information of pairs
  print('Gathering information for pairs')
  clono.count$cdr3 <- GetChainInfo(clono.count$id, 'cdr3')
  clono.count$v_gene <- GetChainInfo(clono.count$id, 'v_gene')
  clono.count$j_gene <- GetChainInfo(clono.count$id, 'j_gene')
  clono.count$chain <- GetChainInfo(clono.count$id, 'chain')
  clono.count$antigen <- GetChainInfo(clono.count$id, 'antigen')
  clono.count$mhc <- GetChainInfo(clono.count$id, 'mhc')
  clono.count$name <- GetPairName(clono.count$id)

  ### Data cleaning
  chain.count$id <- NULL
  chain.count$count <- NULL
  clono.count$id <- NULL
  clono.count$count <- NULL

  ### Formating
  clono.count$index <- lapply(clono.count$index, function(x) as.list(x))
  chain.count$index <- lapply(chain.count$index, function(x) as.list(x))

  return(list(clonoCount = as.list(clono.count),
      chainCount = as.list(chain.count)))
}

GetUnit <- function(study_id, file_logger=NULL) {
  if (is.null(file_logger)) {
    file_logger = log4r::logger()
  }
  log4r::info(file_logger, paste("Getting count unit..."))
  info <- CombineParam(arg$default_params, arg$all_study_params[[study_id]])
  if (!is.null(info$unit) && info$unit != 'unknown') {
    log4r::info(file_logger, paste("Using count unit provided by params..."))
    return(info$unit)
  }
  
  log4r::warn(file_logger, paste0("Count unit provided by params is either null or unknown (", info$unit, "). Attempting to guess unit from count data..."))
  filepath <- GetPath(study_id)
  count_data <- ReadMtxFromBioTuringRawH5(filepath)
  if (sum(abs(round(count_data@x[1:1000], 0) - count_data@x[1:1000])) > 0) { # Non integer
    log4r::warn(file_logger, paste("Count unit is non integer."))
    return("not umi")
  }
  log4r::info(file_logger, paste("Count unit is integer. Assuming that they are UMI"))
  return("umi")  
}

CombineParam <- function(default_params, study_params) {
  # USAGE: Enable user to specify params to be used for a particular study, if differs from the global params
  # Example
  # global_param <- list(dims = 2, perplexity = 10, output.dir = output.dir, seed = 2409, correct_method='harmony')
  # study_param <- list(dims = 3, seed = 2409, correct_method='none')
  # CombineParam(global_param, study_param) --> list(dims = 3, seed = 2409, correct_method='none')
  
  if (is.null(study_params)) {
    study_params = list()
  }
  for (key in names(default_params)) {
    if (is.null(study_params[[key]])) {
      study_params[[key]] <- default_params[[key]]
    }
  }
  return(study_params)
}

RunPipeline <- function(study_id, arg, add_meta = FALSE) {
  RAW_PATH <- arg$raw_path
  OUT_DIR <- arg$output_dir
  OLD_DATA <- arg$old_data

  filepath <- ConnectPath(RAW_PATH, paste0(study_id, '.hdf5'))
  output_path <- ConnectPath(OUT_DIR, paste0(study_id, '.bcs'))
  # default_logger <- log4r::logger()
  log_file <- ConnectPath(OUT_DIR, paste0(study_id, ".log"))
  file_logger <- log4r::logger(appenders = c(log4r::console_appender(), log4r::file_appender(log_file)))
  
  if (file.exists(output_path)) {
    print(paste("This study has already been process:", study_id))

    if (add_meta) {
      log4r::info(file_logger, paste("Adding metadata from file"))
      AddMetadataFromFile(arg, output_path, study_id)
    }
    return(FALSE)
  }

  if (file.exists(log_file)) {
    print(paste("WARNING: log file exist, this might cause some conflict in log: ", log_file))
    log4r::info(file_logger, paste("---------***---------"))
  }
  
  # print(paste("Processing:", study_id))
  log4r::info(file_logger, paste("Processing:", study_id))
  SanityCheck(filepath, file_logger)

  params <- CombineParam(arg$default_params, arg$all_study_params[[study_id]])

  count_data <- ReadMtxFromBioTuringRawH5(filepath)
  
  TrimGeneName <- function(name) {
    if (!IsEnsemblID(name)) { # Only trim for EnsemblID
      return(name)
    }
    if (any(grepl("\\.\\d+$", name))) {
      log4r::warn(file_logger, paste("Trimming suffix (eg `.1`) from EnsemblID..."))
    }
    log4r::info(file_logger, paste("EnsemblID after trimming:"))
    name <- gsub("\\.\\d+$", "", name)
    log4r::info(file_logger, do.call(paste, as.list(name[1:5])))
    return(name)
  }

  ft <- rownames(count_data)
  if (IsEnsemblID(ft)) {
    ft <- TrimGeneName(ft)
    rownames(count_data) <- ft
  }

  if (any(duplicated(colnames(count_data)))) {
    log4r::error(file_logger, paste("There is duplicated barcodes, stop processing."))
    return(FALSE)
  }

  obj <- CreateSeuratObj(count_data, study_id, file_logger)
  
  # TODO: Check this
  unit <- GetUnit(study_id, file_logger)
  if ("ADT" %in% names(obj) && unit == "umi") {
    # print("Normalizing ADT data...")
    log4r::info(file_logger, paste("Normalizing ADT data..."))
    obj <- NormalizeADT(obj)
  }

  # TODO: Handle Clonotype info
  # clonotype_data <- GetClonotypeData(study_id)
  # clonotype_data <- FilterClonotype(clonotype_data, colnames(count_data))
  # CountFromTCRData(clonotype_data)
  meta_data <- GetMetadata(study_id)
  
  if (!is.null(meta_data)) {

    for (meta_name in names(meta_data)) {
      log4r::info(file_logger, paste("Adding metadata:", meta_name))
      meta <- as.numeric(as.character(meta_data[[meta_name]]))
      
      if (any(is.na(meta))) { # To avoid mistakenly forcing character into numeric
        meta <- meta_data[[meta_name]]
      }
      obj <- Seurat::AddMetaData(obj, meta, col.name = meta_name)
    }
  }

  # TODO: Check unit before normalizing?
  # print(paste("Count unit:", unit))
  log4r::info(file_logger, paste("Count Unit", unit))
  if (unit != 'lognorm') {
    log4r::info(file_logger, paste("Normalizing data"))
    obj <- Seurat::NormalizeData(obj)
  } else {
    log4r::info(file_logger, paste("Skip normalization."))
  }

  obj <- Seurat::FindVariableFeatures(obj, nfeatures=min(params$n_variable_features, nrow(obj)))
  obj <- Seurat::ScaleData(obj)
  obj <- Seurat::RunPCA(obj)
  
  correct_method <- GetBatchCorrection(study_id, file_logger)
  key <- 'pca'
  if (correct_method %in% c("cca", "mnn", "harmony")) {
    # print(paste("Run batch correction for:", study_id))
    log4r::info(file_logger, paste("Run HARMONY batch correction for:", study_id))
    batch <- ReadBioTuringRawH5Slot(filepath, "/batch")
    obj <- Seurat::AddMetaData(obj, batch, col.name = 'batch')
    obj <- harmony::RunHarmony(obj, "batch", plot_convergence = FALSE, verbose = TRUE, epsilon.harmony = -Inf, max.iter.harmony = 30)
    key <- 'harmony'
  } else {
    # print(paste("Not correct batch for:", study_id))
    log4r::info(file_logger, paste("Not correct batch for:", study_id))
  }
  
  mat <- obj@reductions[[key]]@cell.embeddings
  mat <- mat[, 1:min(ncol(mat), 50)]
  
  # Run TSNE and UMAP
  log4r::info(file_logger, paste("Running T-SNE and UMAP"))
  tsne.embeddings <- RunTSNE(mat, params)
  obj@reductions[['tSNE']] <- Seurat::CreateDimReducObject(embeddings=tsne.embeddings, assay="RNA", key=paste0(key, "_"))
  
  umap.embeddings <- RunUMAP(mat, params) 
  obj@reductions[['UMAP']] <- Seurat::CreateDimReducObject(embeddings=umap.embeddings, assay="RNA", key=paste0(key, "_"))
  
  # Run Graph-based cluster
  log4r::info(file_logger, paste("Running Graph-based clustering"))
  
  if (params$correct_method %in% c("cca", "mnn", "harmony") && params$correct_method %in% names(obj@reductions)) {
    n.dim.use <- min(30, ncol(obj@reductions[[params$correct_method]]))
    obj <- Seurat::FindNeighbors(obj, dims = 1:n.dim.use, reduction = params$correct_method)
  } else {
    n.dim.use <- min(30, ncol(obj@reductions$pca))
    obj <- Seurat::FindNeighbors(obj, dims = 1:n.dim.use, reduction = "pca")
  }
  
  log4r::info(file_logger, paste("Running Louvain"))
  obj <- Seurat::FindClusters(obj, verbose = FALSE)
  obj@meta.data$bioturing_graph <- as.numeric(obj@meta.data$seurat_clusters)
  obj@meta.data$bioturing_graph <- factor(
    paste("Cluster", obj@meta.data$bioturing_graph),
    levels = paste("Cluster", sort(unique(obj@meta.data$bioturing_graph)))
  )
  
  rhdf5::h5closeAll()
  log4r::info(file_logger, paste("Exporting BCS"))
  rBCS::ExportSeurat(obj, output_path, overwrite=FALSE)
  
  log4r::info(file_logger, paste("Adding metadata from file"))
  AddMetadataFromFile(arg, output_path, study_id)
  
  log4r::info(file_logger, paste("FINISHED"))
}

ProcessBatch <- function(study_list, arg) {
  study_list <- study_list
  lapply(study_list, RunPipeline, arg=arg)
}

GenerateImportArg <- function(study_id, bcs_dir = output_dir, write_dir = "/Users/bioturing/.BioTBData/Data/SingleCell/Study/") {
  
  info <- GetInfo(study_id)
  import_arg <- list()
  import_arg$group_name_aws <- NULL
  import_arg$input_path <- ConnectPath(bcs_dir, paste0(study_id, '.bcs'))
  import_arg$path_type <- "server" # `local` or `server`
  import_arg$batch_names <- NULL
  import_arg$input_id <- c(uuid::UUIDgenerate())
  
  correct_method <- GetBatchCorrection(study_id)
  import_arg$correct_method <- ifelse(correct_method != "none", "harmony", "none") #Study in production could have "cca" "mnn" as correct_method, but for this pipeline, we only use Harmony
  
  import_arg$type <- c("bcs")
  import_arg$unit <- GetUnit(study_id) # TODO: Check this again
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
  # import_arg$db_path <- "_$_DB_MODEL_PATH_$_/"
  import_arg$output_path <- ConnectPath(write_dir, study_id, "main")
  # import_arg$output_path <- ConnectPath("_$_PRIVATE_USER_STUDIES_PATH_$_/vu@bioturing.com", study_id, "main")
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



# TODO: Catch more messages to log.
