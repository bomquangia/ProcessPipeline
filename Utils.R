
TrimGeneName <- function(name, file_logger=NULL) {
  if (!IsEnsemblID(name)) { # Only trim for EnsemblID
    return(name)
  }
  
  if (is.null(file_logger)) {
  file_logger = log4r::logger()
  }
  
  if (any(grepl("\\.\\d+$", name))) {
    log4r::warn(file_logger, paste("Trimming suffix (eg `.1`) from EnsemblID..."))
  }
  log4r::info(file_logger, paste("EnsemblID after trimming:"))
  name <- gsub("\\.\\d+$", "", name)
  log4r::info(file_logger, do.call(paste, as.list(name[1:5])))
  return(name)
}

GetUnit <- function(study_id, file_logger=NULL) {
  if (is.null(file_logger)) {
    file_logger = log4r::logger()
  }
  log4r::info(file_logger, paste("Getting count unit..."))
  
  info <- GetParams(study_id, arg)
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

EstimatePerplexity <- function(n) {
  return(min(max(round(n / 100), 2), 30))
}

GetParams <- function(study_id, arg, harmony_only = TRUE, estimate_perplexity = FALSE, n_data = NULL) {
  params <- CombineParam(arg$default_params, arg$old_study_params[[study_id]])
  
  # All studies requiring batch correction will only use harmony to save computation
  # unless it is specifically overwritten in arg$manual_params
  if (harmony_only) {
    params$correct_method <- ifelse(params$correct_method != "none", "harmony", "none")
  }
  
  if (estimate_perplexity) {
    print("Estimating perlexity...")
    print(paste("Before:", params$perplexity))
    params$perplexity <- EstimatePerplexity(n_data)
    print(paste("After:", params$perplexity))
  }

  params <- CombineParam(params, arg$manual_params[[study_id]])
  return(params)
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

ConnectPath <- function(...) {
  res <- file.path(...)
  if (.Platform$OS.type == 'windows') {
    res <- gsub('\\/', '\\\\', res)
  }
  return(res)
}

SanityCheck <- function(filepath, file_logger=NULL) {
  if (is.null(file_logger)) {
    file_logger = log4r::logger()
  }
  dims = ReadBioTuringRawH5Slot(filepath, "/shape")
  
  features <- ReadBioTuringRawH5Slot(filepath, "/features")
  barcodes <- ReadBioTuringRawH5Slot(filepath, "/barcodes")

  n <- length(features)
  print(features[c(1:5, (n-5):n)])
  
  log4r::info(file_logger, paste("Example features: "))
  log4r::info(file_logger, paste(do.call(paste, as.list(features[c(1:5, (n-5):n)]))))
  
  log4r::info(file_logger, paste("Example barcodes: "))
  log4r::info(file_logger, paste(do.call(paste, as.list(barcodes[1:5]))))

  batch <- ReadBioTuringRawH5Slot(filepath, "/batch")

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

HexToRaw <- function(text) {
  vals <- matrix(as.integer(as.hexmode(strsplit(text, "")[[1]])), ncol=2, byrow=TRUE)
  vals <- vals %*% c(16, 1)
  return(as.raw(vals))
}

WriteJSON <- function(object, filepath, encrypt=FALSE, ...) {
  CreateDir(dirname(filepath))
  if (nchar(Sys.getenv('e32dc2')) > 0 && encrypt) {
    EncryptJSON(object, filepath, ...)
  } else {
    jsonlite::write_json(object, filepath, ...)
  }
}

EncryptJSON <- function(object, filepath, ...) {
  secret.key <- Sys.getenv('e32dc2')
  iv <- HexToRaw(digest::digest(runif(1), algo="md5"))
  aes <- digest::AES(charToRaw(secret.key), mode="CTR", IV=iv)
  raw.content <- charToRaw(jsonlite::toJSON(object, ...))
  content <- suppressWarnings(aes$encrypt(raw.content))
  result <- list(iv = paste(iv, collapse=""), content = paste(content, collapse=""))
  jsonlite::write_json(result, path=filepath, auto_unbox=TRUE)
}


IsBioTuring <- function(x) {
  if (is.null(x)) {
    return(FALSE)
  }
  return(grepl("@bioturing[.]com$", x))
}


ReadJSON <- function(filepath, ...) {
  
  secret.key <- Sys.getenv('e32dc2')
  code <- try(jsonlite::fromJSON(filepath))
  
  if (inherits(code, "try-error")) {
    print("cannot read run_info.json from zip, consider do it manually")
    return(NULL)
  }
  
  if (!'iv' %in% names(code)) { # For unencrypted run_info.json
    return(code)
  }
  aes <- digest::AES(charToRaw(secret.key), mode="CTR", IV=HexToRaw(code$iv))
  raw.content <- suppressWarnings(aes$decrypt(HexToRaw(code$content)))
  return(jsonlite::fromJSON(raw.content, ...))
}

GenerateImportArg <- function(study_id, bcs_dir = output_dir, write_dir = "/Users/bioturing/.BioTBData/Data/SingleCell/Study/") {
  
  info <- GetInfo(study_id)
  import_arg <- list()
  import_arg$group_name_aws <- NULL
  import_arg$input_path <- ConnectPath(bcs_dir, paste0(study_id, '.bcs'))
  import_arg$path_type <- "server" # `local` or `server`
  import_arg$batch_names <- NULL
  import_arg$input_id <- c(uuid::UUIDgenerate())
  
  import_arg$correct_method <- GetBatchCorrection(study_id)
  
  import_arg$type <- c("bcs")
  import_arg$unit <- GetUnit(study_id) # TODO: Check this again
  import_arg$email <- "vu@bioturing.com"
  import_arg$species <- info$species
  import_arg$log_path <- ConnectPath(write_dir, study_id, "tmp/submit.log")
  import_arg$dimred_method <- c("tsne", "umap")
  import_arg$dimred_perplexity <- 0 #fake
  import_arg$dimred_mode <- "slow" #fake
  import_arg$quant <- list(method = "unknown", ref = "unknown")
  import_arg$seed <- 2409
  import_arg$platform <- "unknown"
  import_arg$filter <- list(cell=0, gene = c(0,0), mito = 100, top = 2000) #fake
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

RunDiagnostics <- function(study_id, arg, file_logger=NULL) {
  if (is.null(file_logger)) {
    file_logger = log4r::logger()
  }
  
  log4r::info(file_logger, paste("Running Diagnosis for", study_id))
  study_id <- as.character(study_id)
  bcs_path <- ConnectPath(arg$output_dir, paste0(study_id, '.bcs'))
  hdf5_path <- ConnectPath(arg$raw_path, paste0(study_id, '.hdf5'))
  
  old_study_path <- ConnectPath(arg$old_data, study_id)

  # Info needed to collect: Processed, Same barcodes, Number of batches, batch_correct, 
  processed <- file.exists(bcs_path)
  
  hdf5_exists <- file.exists(hdf5_path)
  
  old_study_exists <- dir.exists(old_study_path)
  
  barcodes <- try(ReadBioTuringRawH5Slot(hdf5_path, "/barcodes"))
  if (inherits(barcodes, "try-error")) {
    log4r::info(file_logger, paste("Cannot open hdf5"))
    return(NULL)
  }
  
  studyInfo <- GetInfo(study_id)

  original_n_batch <- length(unique(studyInfo$batch))
  hasADT = sum(studyInfo$ADT_indices) > 0
  filter <- list(gene = c(0,0), mito=100, top =2000, cell=0)
  if (!old_study_exists) {
    log4r::info(file_logger, paste("Either Old study zip not available, or unzip went wrong"))

    return(list(
      study_id=study_id,
      processed=processed,
      hdf5_exists=hdf5_exists,
      original_n_batch=original_n_batch,
      old_n_batch= NULL,
      correct_method= NULL,
      norm_method = NULL,
      hasADT = hasADT,
      unit = NULL,
      filter = NULL))
  }
  
  run_info <- ReadJSON(ConnectPath(arg$old_data, study_id, "run_info.json"))
  
  old_n_batch <- run_info$n_batch
  correct_method <- run_info$ana_setting$batchRemoval
  norm_method <- run_info$ana_setting$normMethod
  unit <- run_info$unit
  
  if (is.null(unit)) {
    unit <- run_info$ana_setting$unit
  }
  
  filter <- GetAnalysisSettings(study_id, arg)$filter

  return(list(
    study_id=study_id,
    processed = processed,
    hdf5_exists=hdf5_exists,
    original_n_batch=original_n_batch,
    old_n_batch=old_n_batch,
    correct_method=correct_method,
    norm_method = norm_method,
    hasADT = hasADT,
    unit = unit,
    filter = filter)
    )
}

GetAnalysisSettings <- function(study_id, arg) {
  # Get run_info$ana_settings from old studies
  old_study_path <- ConnectPath(arg$old_data, study_id)
  if (!dir.exists(old_study_path)) {
    print("Old study folder not found, please make sure that it is present at `arg$old_data`/study_id")
    return(NULL)
  }
  run_info <- ReadJSON(ConnectPath(old_study_path, 'run_info.json'))
  return(run_info$ana_setting)
}

DiagnoseBatch <- function(output_dir, raw_path, old_data, study_list) {
  old_data <- old_data
  output_dir <- output_dir
  raw_path <- raw_path
  diagnosis <- lapply(study_list, RunDiagnostics)

  return(do.call(rbind, diagnosis))
}

CollectOldParams <- function(study_id, arg) {
  info <- RunDiagnostics(study_id, arg = arg)
  ana_settings <- GetAnalysisSettings()
  return(list(correct_method = info$correct_method, norm_method = info$norm_method, unit=info$unit, old_n_batch=info$old_n_batch))
}

GetPrefix <- function(barcode) {
  locations <- stringr::str_locate(barcode, "[ACGT]{16}")
  prefix <- substr(barcode, 1, locations[1]-1)
  return(prefix)
}

GetSuffix <- function(barcode) {
  barcode <- as.character(barcode)
  locations <- stringr::str_locate(barcode, "[ACGT]{16}")
  suffix <- substr(barcode, locations[2]+1, nchar(barcode))
  return(suffix)
}

GetBarcodeCore <- function(barcode) {
  locations <- stringr::str_locate(barcode, "[ACGT]{16}")
  core <- substr(barcode, locations[1], locations[2])
  return(core)
}

CheckUnique <- function(barcodes) {
  return(length(unique(barcodes)) == length(barcodes))
}

CheckBarcodes <- function(study_id) {
  print(paste0("Checking study: ", study_id))
  new_study <- ConnectPath(output_dir, study_id)
  old_study <- ConnectPath(old_data, study_id) # Should change this to read from 

  if (!file.exists(new_study)) {
    print('Not yet finished or you are pointing to the wrong directory for new study')
    return(c(NA,NA,NA,NA))
  }
  if (!file.exists(old_study)) {
    print('No data of old study available or you are pointing to the wrong directory for new study')
    return(c(study_id,NA,NA,NA))
  }

  new_barcodes <- as.character(read.csv(ConnectPath(new_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
  old_barcodes <- as.character(read.csv(ConnectPath(old_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
  
  match <- old_barcodes %in% new_barcodes
  full_match <- paste("Full match: ", sum(match), '/', length(match))
  print(full_match)
  return(list(study_id = study_id, Number_Match= sum(match), Old_barcodes = length(old_barcodes), New_barcodes = length(new_barcodes)))
}



