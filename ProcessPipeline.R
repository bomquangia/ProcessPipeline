source("/mnt2/vu/script/ProcessPipeline/AddClonotype.R")
source("/mnt2/vu/script/ProcessPipeline/BatchCorrection.R")
source("/mnt2/vu/script/ProcessPipeline/FindClusters.R")
source("/mnt2/vu/script/ProcessPipeline/ImportMetadata.R")
# source("/mnt2/vu/script/ProcessPipeline/MapOldMetaData.R")
source("/mnt2/vu/script/ProcessPipeline/MatchBarcodes.R")
source("/mnt2/vu/script/ProcessPipeline/ReduceDimension.R")
source("/mnt2/vu/script/ProcessPipeline/Utils.R")


RunPipeline <- function(study_id, arg, add_meta = FALSE) {
  RAW_PATH <- arg$raw_path
  OUT_DIR <- arg$output_dir
  OLD_DATA <- arg$old_data

  filepath <- ConnectPath(RAW_PATH, paste0(study_id, '.hdf5'))
  output_path <- ConnectPath(OUT_DIR, paste0(study_id, '.bcs'))

  log_file <- ConnectPath(OUT_DIR, paste0(study_id, ".log"))
  file_logger <- log4r::logger(appenders = c(log4r::console_appender(), log4r::file_appender(log_file)))
  
  if (file.exists(output_path)) {
    print(paste("This study has already been process:", study_id))

    if (add_meta) {
      log4r::info(file_logger, paste("Adding metadata from file"))
      AddMetadataFromFile(arg, output_path, study_id, file_logger)
    }
    return(FALSE)
  }

  if (file.exists(log_file)) {
    print(paste("WARNING: log file exist, this might cause some conflict in log: ", log_file))
    log4r::info(file_logger, paste("---------***---------"))
  }
  
  log4r::info(file_logger, paste("Processing:", study_id))
  SanityCheck(filepath, file_logger)

  params <- GetParams(study_id, arg)

  count_data <- ReadMtxFromBioTuringRawH5(filepath)
  
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
  
  params$unit <- GetUnit(study_id, file_logger)
  log4r::info(file_logger, paste("Count Unit:", params$unit))
  # 
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

  params$correct_method <- GetBatchCorrection(study_id, file_logger)

  batch <- ReadBioTuringRawH5Slot(filepath, "/batch")
  obj <- Seurat::AddMetaData(obj, batch, col.name = 'bioturing_batch')

  obj <- Preprocess(obj, params, file_logger, arg$seed)
  obj <- ReduceDimension(obj, params, file_logger, arg$seed)
  obj <- FindClusters(obj, file_logger, params)
  # obj@meta.data <- data.frame(
  #                   bioturing_graph = paste("Cluster", as.numeric(obj@meta.data$seurat_clusters)))
  
  rhdf5::h5closeAll()
  log4r::info(file_logger, paste("Exporting BCS"))
  rBCS::ExportSeurat(obj, output_path, overwrite=FALSE)
  
  log4r::info(file_logger, paste("Adding metadata from file"))
  AddMetadataFromFile(arg, output_path, study_id, file_logger)
  
  log4r::info(file_logger, paste("Generating import argument for NoraSC::ProcessFromObject"))
  import_arg <- GenerateImportArg(study_id, bcs_dir = OUT_DIR)
  jsonlite::write_json(import_arg, ConnectPath(OUT_DIR, paste0(study_id, '.json')))

  log4r::info(file_logger, paste("FINISHED"))
}

ProcessBatch <- function(study_list, arg) {
  study_list <- study_list
  lapply(study_list, RunPipeline, arg=arg)
}
