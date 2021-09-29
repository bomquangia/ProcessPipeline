ImportMetadata <- function(arg.json, file_logger) {
  
  IsNumeric <- function(x) {
    y <- as.numeric(as.character(x))
    return((!any(is.na(y))) && length(unique(x)) >= 100)
  }

  MapBarcode <- function(barcodes, meta) {
    # Filter cells that are not in the dataset
    index.matches <- MatchBarcodes(meta[[1]], barcodes)
    
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
    special.names <- c("Graph-based clusters", "Graph based", "Batch", "Total count", "Total expressed feature")
    data <- ReadPlainTextTable(filepath)
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
    
    omit_columns <- grepl("Barcodes", colnames(meta), ignore.case = TRUE)
    omit_columns[1] <- FALSE
    
    meta <- meta[, !omit_columns]
    meta <- MapBarcode(barcodes, meta)
    log4r::info(file_logger, paste("Name of metadata to be import: "))
    log4r::info(file_logger, do.call(paste, as.list(names(meta))))
    
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

AddMetadataFromFile <- function(arg, output_path, study_id, file_logger=NULL) {
  if (is.null(file_logger)) {
    file_logger = log4r::logger()
  }
  setwd(arg$output_dir)
  unzip(output_path, exdir = arg$output_dir, unzip="unzip")
  old_dir <- unzip(output_path, list=TRUE)$Name[1]
  import_arg <- list(input_path=ConnectPath(arg$meta_dir, paste0(study_id, '.tsv')),
                  data_path=ConnectPath(old_dir, "/main"),
                  # raw_h5 = GetPath(study_id),
                  output_path= ConnectPath(arg$output_dir, old_dir, 'main/metadata'),
                  type = "single",
                  email = "vu@bioturing.com",
                  unique_limit = 100,
                  bbrowser_version = "2.9.23")
  ImportMetadata(jsonlite::toJSON(import_arg), file_logger)
  zip::zip(paste0(study_id, '.bcs'), old_dir, compression_level=1)
  unlink(old_dir, recursive = TRUE)
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

GenerateUUID <- function() {
  return (gsub('-', '', uuid::UUIDgenerate()))
}

CreateHistory <- function(options) {
  if (is.null(options$hash_id)) {
    options$hash_id <- GenerateUUID()
  }
  if (is.null(options$description)) {
    options$description <- "Created by BBrowser"
  }
  content <- list(list(
    created_by = options$email,
    created_at = as.numeric(Sys.time()) * 1000,
    hash_id = options$hash_id,
    description = options$description
  ))
  return(content)
}

DetectDelim <- function(filepath) {
  if (endsWith(filepath, ".tsv")) {
    return("\t")
  } else if (endsWith(filepath, ".csv")) {
    return(",")
  }

  lines <- readLines(filepath, n = 50)
  delim.list <- c("\t", ",", " ")
  len.list <- lapply(delim.list, function(delim) {
    return(sapply(strsplit(lines, delim), length))
  })
  is.valid <- sapply(len.list, function(x) {
    if (all(x[2:length(x)] == x[2])) {
      dif <- abs(x[1] - x[2])
      if (dif <= 1) {
        return(x[2])
      }
    }
    return(-1)
  })
  if (is.valid[1] > 1) { # Favor tabs
    is.valid[-1] <- -1
  }
  return(delim.list[order(is.valid, decreasing = TRUE)][1])
}

ReadPlainTextTable <- function(tb.path, delim=NULL, convert=TRUE, ...) {
  cols <- readr::cols(.default = readr::col_character()) # Force everything to strings
  na.chars <- c("", character(0)) # There will be recognize as NA
  delim <- if (is.null(delim)) DetectDelim(tb.path) else delim
  data <- readr::read_delim(tb.path, delim=delim, na=na.chars, col_types=cols, ...)
  if (convert) {
    data <- readr::type_convert(data) # auto detect column type
  }
  return(data)
}

CreateDir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

