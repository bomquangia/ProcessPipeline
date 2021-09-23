

## GSE123814_BCC
study_id <- "GSE123814_BCC"
output_dir <- "/mnt2/vu/script/output"
old_data <- "/mnt2/vu/script/old_data"
new_study <- ConnectPath(output_dir, study_id)
old_study <- ConnectPath(old_data, study_id) # Should change this to read from 

new_barcodes <- as.character(read.csv(ConnectPath(new_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
old_barcodes <- as.character(read.csv(ConnectPath(old_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)

new_barcodes[1:5]
old_barcodes[1:5]

# --> Two set of barcodes mismatch by . and _ --> Need to substitute all `.` in old_barcodes into `_` to get new_barcode
# First, need to check that no barcode in new_barcodes contain `.`
any(grepl('.', new_barcodes, fixed=TRUE)) # --> FALSE --> Good to go

# Substuting dash for dot
modified_old_barcodes <- gsub('.','_', old_barcodes, fixed=TRUE)

# Check that everything fit perfectly
all(new_barcodes %in% modified_old_barcodes) # --> TRUE --> Good to go

# Generating tsv file containing old_barcodes new_barcodes and all old metadata
old_meta_file <- read.csv("/mnt2/vu/script/metadata/GSE123814_BCC/GSE123814_BCC_metadata_1632378535407.tsv", sep="\t")
barcodes_mapping <- data.frame(cbind(old_barcodes, new_barcodes))

# Merging by colnames
old_meta_file <- 
merged_meta <- merge(old_meta_file, barcodes_mapping, by.x=c("Barcodes"), by.y=c("old_barcodes"))
write.table(merged_meta, "/mnt2/vu/script/metadata/GSE123814_BCC/GSE123814_BCC_merged_meta.tsv", row.names = FALSE, sep="\t")

## ------- ## 

## GSE123814_SCC

study_id <- "GSE123814_SCC"
output_dir <- "/mnt2/vu/script/output"
old_data <- "/mnt2/vu/script/old_data"
new_study <- ConnectPath(output_dir, study_id)
old_study <- ConnectPath(old_data, study_id) # Should change this to read from 

new_barcodes <- as.character(read.csv(ConnectPath(new_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
old_barcodes <- as.character(read.csv(ConnectPath(old_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)

new_barcodes[1:5]
old_barcodes[1:5]

# --> Two set of barcodes mismatch by . and _ --> Need to substitute all `.` in old_barcodes into `_` to get new_barcode
# First, need to check that no barcode in new_barcodes contain `.`
stopifnot(!any(grepl('.', new_barcodes, fixed=TRUE))) # --> Good to go

# Substuting dash for dot
modified_old_barcodes <- gsub('.','_', old_barcodes, fixed=TRUE)

# Check that everything fit perfectly
stopifnot(all(new_barcodes %in% modified_old_barcodes)) # --> TRUE --> Good to go

# Generating tsv file containing old_barcodes new_barcodes and all old metadata
old_meta_file <- read.csv("/mnt2/vu/script/metadata/GSE123814_SCC/GSE123814_SCC_metadata_1632379587751.tsv", sep="\t")
barcodes_mapping <- data.frame(cbind(old_barcodes, new_barcodes))

# Merging by colnames
merged_meta <- merge(old_meta_file, barcodes_mapping, by.x=c("Barcodes"), by.y=c("old_barcodes"))
write.table(merged_meta, "/mnt2/vu/script/metadata/GSE123814_SCC/GSE123814_SCC_merged_meta.tsv", row.names = FALSE, sep="\t")


# GSE123904_human
study_id <- "GSE123904_human"
output_dir <- "/mnt2/vu/script/output"
old_data <- "/mnt2/vu/script/old_data"
jnj_meta <- "/mnt1/tracy/jnj/metadata"
metadata_dir <- ConnectPath("/mnt2/vu/script/metadata", study_id)

new_study <- ConnectPath(output_dir, study_id)
old_study <- ConnectPath(old_data, study_id) # Should change this to read from 

if (!dir.exists(metadata_dir)) {
  dir.create(metadata_dir)
}

meta_path <- ConnectPath(jnj_meta, paste0(study_id, '.tsv'))
if (file.exists(meta_path)) { # JnJ Metadata exists, take old barcodes from there
  old_barcodes <- read.csv(meta_path, sep ="\t")$Barcodes
} else { # Export metadata manually and read it
  old_barcodes <- as.character(read.csv(ConnectPath(old_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
}

new_barcodes <- as.character(read.csv(ConnectPath(new_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)

new_barcodes[1:5]
old_barcodes[1:5]
# "MSK_LX255B_METASTASIS_120703408781739" vs "1" --> Ummappable



# Generating tsv file containing old_barcodes new_barcodes and all old metadata
old_meta_file <- read.csv("/mnt2/vu/script/metadata/GSE123814_SCC/GSE123814_SCC_metadata_1632379587751.tsv", sep="\t")
# barcodes_mapping <- data.frame(cbind(old_barcodes, new_barcodes))