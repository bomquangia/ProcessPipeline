CheckBarcodes <- function(new.study, old.study) {
  # new.study and old.study are path to respective folder
  # browser()
  new.barcodes <- read.csv(ConnectPath(new.study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1
  old.barcodes <- read.csv(ConnectPath(old.study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1
  if (length(new.barcodes) != length(old.barcodes)) {
    print('Different number of cells')
    return(FALSE)
  }
  
  return(all(as.character(new.barcodes)==as.character(old.barcodes)))
}

# Test code
study <- "GSE123022"
old.data <- "/mnt2/vu/script/old_data"
new.data <- "/mnt2/vu/script/output"

CheckBarcodeWrapper <- function(study) {
  print(paste0("Checking study: ", study))
  # if (study == "GSE123814_BCC") {
  #   browser()
  # }
  new.study <- ConnectPath(new.data, study)
  old.study <- ConnectPath(old.data, study)
  if (!file.exists(new.study)) {
    return('Not yet finished')
  }
  CheckBarcodes(new.study, old.study)
}
match <- lapply(study_list$V1, CheckBarcodeWrapper)