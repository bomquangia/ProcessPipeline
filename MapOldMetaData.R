
CheckBarcodes <- function(new_study, old_study) {
  # new.study and old.study are path to respective folder
  new_barcodes <- as.character(read.csv(ConnectPath(new_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
  old_barcodes <- as.character(read.csv(ConnectPath(old_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
  
  match <- old_barcodes %in% new_barcodes
  full_match <- paste("Full match: ", sum(match), '/', length(match))
  
  return(full_match)
}

CheckBarcodeWrapper <- function(study) {
  print(paste0("Checking study: ", study))

  new_study <- ConnectPath(output_dir, study)
  old_study <- ConnectPath(old_data, study) # Should change this to read from 
  if (!file.exists(new_study)) {
    return('Not yet finished')
  }
  CheckBarcodes(new_study, old_study)
}

HexToRaw <- function(text) {
  vals <- matrix(as.integer(as.hexmode(strsplit(text, "")[[1]])), ncol=2, byrow=TRUE)
  vals <- vals %*% c(16, 1)
  return(as.raw(vals))
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

RunDiagnostics <- function(study_id, file_logger=NULL) {
  if (is.null(file_logger)) {
    file_logger = log4r::logger()
  }
  
  log4r::info(file_logger, paste("Running Diagnosis for", study_id))
  study_id <- as.character(study_id)
  bcs_path <- ConnectPath(output_dir, paste0(study_id, '.bcs'))
  hdf5_path <- ConnectPath(raw_path, paste0(study_id, '.hdf5'))
  
  old_study_path <- ConnectPath(old_data, study_id)

  # Info needed to collect: Processed, Same barcodes, Number of batches, batch_correct, 
  processed <- file.exists(bcs_path)
  
  hdf5_exists <- file.exists(hdf5_path)
  # browser()
  old_study_exists <- dir.exists(old_study_path)
  
  barcodes <- try(readThaoH5Slot(hdf5_path, "/barcodes"))
  if (inherits(barcodes, "try-error")) {
    log4r::info(file_logger, paste("Cannot open hdf5"))
    return(NULL)
  }
  
  studyInfo <- GetInfo(study_id)

  original_n_batch <- length(unique(studyInfo$batch))
  hasADT = sum(studyInfo$ADT_indices) > 0
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
      unit = NULL))
  }
  # same_barcodes <- CheckBarcodeWrapper(study_id)
  
  run_info <- ReadJSON(ConnectPath(old_data, study_id, "run_info.json"))
  
  old_n_batch <- run_info$n_batch
  correct_method <- run_info$ana_setting$batchRemoval
  norm_method <- run_info$ana_setting$normMethod
  unit <- run_info$unit
  
  if (is.null(unit)) {
    unit <- run_info$ana_setting$unit
  }

  return(list(
    study_id=study_id,
    processed = processed,
    hdf5_exists=hdf5_exists,
    original_n_batch=original_n_batch,
    old_n_batch=old_n_batch,
    correct_method=correct_method,
    norm_method = norm_method,
    hasADT = hasADT,
    unit = unit)
    )
}

DiagnoseBatch <- function(output_dir, raw_path, old_data, study_list) {
  old_data <- old_data
  output_dir <- output_dir
  raw_path <- raw_path
  diagnosis <- lapply(study_list, RunDiagnostics)

  return(do.call(rbind, diagnosis))
}


### Run diagnosis 
# output_dir <- '/mnt2/vu/script/output'
# raw_path <- "/mnt2/bioturing_data/raw_hdf5/"
# old_data <- "/mnt2/vu/script/old_data"
# batch_4 = "/mnt2/bioturing_data/files/batch_4.tsv"
# study_list = read.table(file = batch_4, sep = '\t', header = FALSE)$V1

# diagnosis <- DiagnoseBatch(output_dir, raw_path, old_data, study_list)

## Collect Old run_info.json and compile a list of param
CollectOldParams <- function(study_id) {
  info <- RunDiagnostics((study_id))
  # params <- list()
  # params[[study_id]] <- list(correct_method = info$correct_method, norm_method = info$norm_method, unit=info$unit, old_n_batch=info$old_n_batch)
  return(list(correct_method = info$correct_method, norm_method = info$norm_method, unit=info$unit, old_n_batch=info$old_n_batch))
}
