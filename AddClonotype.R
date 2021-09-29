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
