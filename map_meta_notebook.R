

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
jnj_meta <- "/mnt1/tracy/jnj/metadata"
old_meta_file <- read.csv(ConnectPath(jnj_meta, paste0(study_id, ".tsv")), sep="\t")
barcodes_mapping <- data.frame(cbind(old_barcodes, new_barcodes))
# barcodes_mapping <- lapply(barcodes_mapping, as.factor)
# Merging by colnames

merged_meta <- merge(old_meta_file, barcodes_mapping, by.x=c("Barcodes"), by.y=c("old_barcodes"))
merged_meta[1] <- merged_meta$new_barcodes
merged_meta <- cbind(merged_meta$new_barcodes, merged_meta)
colnames(merged_meta)[1] <- "h5_barcodes"
write.table(merged_meta, "/mnt2/vu/script/metadata/GSE123814_BCC/GSE123814_BCC.tsv", row.names = FALSE, sep="\t")

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
old_meta_file <- read.csv(ConnectPath("GSE123814_SCC_metadata_1632379587751.tsv", sep="\t")
barcodes_mapping <- data.frame(cbind(old_barcodes, new_barcodes))

# Merging by colnames
merged_meta <- merge(old_meta_file, barcodes_mapping, by.x=c("Barcodes"), by.y=c("old_barcodes"))
write.tsv(merged_meta, "/mnt2/vu/script/metadata/GSE123814_SCC/GSE123814_SCC_merged_meta.tsv", row.names = FALSE)


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
  print("JnJ metadata exists")
  old_meta_file <- read.csv(meta_path, sep ="\t")
  barcode_slot <- ifelse('BBrowser_barcodes' %in% colnames(old_meta_file), 'BBrowser_barcodes', "Barcodes")
  old_barcodes <- as.character(old_meta_file[[barcode_slot]])
} else { # Export metadata manually and read it
  print("JnJ metadata not exists, take barcodes.tsv from old study zip")
  old_barcodes <- as.character(read.csv(ConnectPath(old_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
}

# production_metadata <- read.csv("/mnt2/vu/script/metadata/GSE124310/GSE124310_metadata_1632386464074.tsv", sep="\t")
new_barcodes <- as.character(read.csv(ConnectPath(new_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)

new_barcodes[1:5] #"H15w_Z1_AAACCTGAGAACTCGG-1"
old_barcodes[1:5] #1_AAACCTGAGGGTGTGT

# "MSK_LX255B_METASTASIS_120703408781739" vs "1" --> Ummappable

# GSE124310
study_id <- "GSE124310"
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

production_metadata <- read.csv("/mnt2/vu/script/metadata/GSE124310/GSE124310_metadata_1632386464074.tsv", sep="\t")
new_barcodes <- as.character(read.csv(ConnectPath(new_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)

new_barcodes[1:5] #"MGUS-1_138N_AAACCTGAGAGGTTGC-1"
old_barcodes[1:5] #AAACGGGCAAGGTTTC

# Need to parse the barcode part and check if unique 

truncated_new_barcodes <- stringr::str_extract(new_barcodes, "[ACGT]{16}")
stopifnot(length(unique(truncated_new_barcodes)) == length(truncated_new_barcodes))
length(unique(old_barcodes))

# GSE124335: 1696/6144 --> No. of cells too different to map
# GSE124395: 24174/60864 --> No. of cells too different to map

# 11. GSE124472_15wk
# 
study_id <- "GSE124472_15wk"
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
  print("JnJ metadata exists")
  old_meta_file <- read.csv(meta_path, sep ="\t")
  barcode_slot <- ifelse('BBrowser_barcodes' %in% colnames(old_meta_file), 'BBrowser_barcodes', "Barcodes")
  old_barcodes <- as.character(old_meta_file[[barcode_slot]])
} else { # Export metadata manually and read it
  print("JnJ metadata not exists, take barcodes.tsv from old study zip")
  old_barcodes <- as.character(read.csv(ConnectPath(old_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
}

# production_metadata <- read.csv("/mnt2/vu/script/metadata/GSE124310/GSE124310_metadata_1632386464074.tsv", sep="\t")
new_barcodes <- as.character(read.csv(ConnectPath(new_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)

new_barcodes[1:5] #"H15w_Z1_AAACCTGAGAACTCGG-1"
old_barcodes[1:5] #1_AAACCTGAGGGTGTGT

# Need to parse the barcode part and check if unique 
truncated_new_barcodes <- stringr::str_extract(new_barcodes, "[ACGT]{16}")
length(unique(truncated_new_barcodes)) == length(truncated_new_barcodes)
stopifnot(length(unique(old_barcodes)) == length(old_barcodes))

# Try to change prefix `H15w_Z1_` --> 1_ and `H15w_Z2_` --> 2
new_prefixes <- sapply(new_barcodes, GetPrefix) # Two unique prefixes
new_suffixes <- sapply(new_barcodes, GetSuffix) # 1 unique suffixes, ignore
new_cores <- sapply(new_barcodes, GetBarcodeCore)

modified_barcodes <- sapply(seq_along(old_barcodes), function(i) {
  prefix <- GetPrefix(old_barcodes[i])
  core <- GetBarcodeCore(old_barcodes[i])
  # suffix <- GetSuffix(old_barcodes[i])
  if (prefix == "1_") {
    modified_prefix <- "H15w_Z1_"
  } else if (prefix == "2_") {
    modified_prefix <- "H15w_Z2_"
  }
  modified_barcode <- paste0(modified_prefix, core, "-1")
  print(paste("Modified barcode:", modified_barcode))
  stopifnot(modified_barcode %in% new_barcodes)
  return(modified_barcode)
})

stopifnot(!any(duplicated(modified_barcodes))) # No duplicate

# With all the above tests passed, we can probably assume that `H15w_Z1_` --> 1_ and `H15w_Z2_` --> 2.
# Now map. 
barcodes_mapping <- cbind(old_barcodes, modified_barcodes)

merged_meta <- merge(old_meta_file, barcodes_mapping, by.x=c(barcode_slot), by.y=c("old_barcodes"))
merged_meta_path <- ConnectPath(metadata_dir, paste0(study_id, "_merged_meta.tsv"))
write.table(merged_meta, merged_meta_path, row.names = FALSE, sep="\t")


# 12. GSE124472_17wk
# 
study_id <- "GSE124472_17wk"
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
  print("JnJ metadata exists")
  old_meta_file <- read.csv(meta_path, sep ="\t")
  barcode_slot <- ifelse('BBrowser_barcodes' %in% colnames(old_meta_file), 'BBrowser_barcodes', "Barcodes")
  old_barcodes <- as.character(old_meta_file[[barcode_slot]])
} else { # Export metadata manually and read it
  print("JnJ metadata not exists, take barcodes.tsv from old study zip")
  old_barcodes <- as.character(read.csv(ConnectPath(old_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
}

# production_metadata <- read.csv("/mnt2/vu/script/metadata/GSE124310/GSE124310_metadata_1632386464074.tsv", sep="\t")
new_barcodes <- as.character(read.csv(ConnectPath(new_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)

new_barcodes[1:5] #"H15w_Z1_AAACCTGAGAACTCGG-1"
old_barcodes[1:5] #1_AAACCTGAGGGTGTGT

# Need to parse the barcode part and check if unique 
truncated_new_barcodes <- stringr::str_extract(new_barcodes, "[ACGT]{16}")
length(unique(truncated_new_barcodes)) == length(truncated_new_barcodes)
stopifnot(length(unique(old_barcodes)) == length(old_barcodes))

# Try to change prefix `H15w_Z1_` --> 1_ and `H15w_Z2_` --> 2
new_prefixes <- sapply(new_barcodes, GetPrefix) # Two unique prefixes
new_suffixes <- sapply(new_barcodes, GetSuffix) # 1 unique suffixes, ignore
new_cores <- sapply(new_barcodes, GetBarcodeCore)

modified_barcodes <- sapply(seq_along(old_barcodes), function(i) {
  prefix <- GetPrefix(old_barcodes[i])
  core <- GetBarcodeCore(old_barcodes[i])
  # suffix <- GetSuffix(old_barcodes[i])
  if (prefix == "1_") {
    modified_prefix <- "H17w_Z1_"
  } else if (prefix == "2_") {
    modified_prefix <- "H17w_Z2_"
  }
  modified_barcode <- paste0(modified_prefix, core, "-1")
  print(paste("Modified barcode:", modified_barcode))
  stopifnot(modified_barcode %in% new_barcodes)
  print("Checking if modified barcode is in new_barcodes")
  # print(modified_barcode %in% new_barcodes)
  return(modified_barcode)
})

stopifnot(!any(duplicated(modified_barcodes))) # No duplicate

# With all the above tests passed, we can probably assume that `H15w_Z1_` --> 1_ and `H15w_Z2_` --> 2.
# Now map. 
barcodes_mapping <- cbind(old_barcodes, modified_barcodes)

merged_meta <- merge(old_meta_file, barcodes_mapping, by.x=c(barcode_slot), by.y=c("old_barcodes"))
merged_meta_path <- ConnectPath(metadata_dir, paste0(study_id, "_merged_meta.tsv"))
write.table(merged_meta, merged_meta_path, row.names = FALSE, sep="\t")



#GSE124472_kidney_organoid

study_id <- "GSE124472_kidney_organoid"
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
  print("JnJ metadata exists")
  old_meta_file <- read.csv(meta_path, sep ="\t")
  barcode_slot <- ifelse('BBrowser_barcodes' %in% colnames(old_meta_file), 'BBrowser_barcodes', "Barcodes")
  old_barcodes <- as.character(old_meta_file[[barcode_slot]])
} else { # Export metadata manually and read it
  print("JnJ metadata not exists, take barcodes.tsv from old study zip")
  old_barcodes <- as.character(read.csv(ConnectPath(old_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
}

# production_metadata <- read.csv("/mnt2/vu/script/metadata/GSE124310/GSE124310_metadata_1632386464074.tsv", sep="\t")
new_barcodes <- as.character(read.csv(ConnectPath(new_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)

new_barcodes[1:5] #"HuOrg_D16_1_AAACCTGAGCGTAATA-1"
old_barcodes[1:5] #1_AAACCTGAGCGTAATA

# Need to parse the barcode part and check if unique 
truncated_new_barcodes <- stringr::str_extract(new_barcodes, "[ACGT]{16}")
length(unique(truncated_new_barcodes)) == length(truncated_new_barcodes)
stopifnot(length(unique(old_barcodes)) == length(old_barcodes))

# Try to change prefix  "HuOrg_D16_1_" --> "_1" and  "HuOrg_D16_2_" --> "_2" and "HuOrg_D28_1_" --> "_3" and "HuOrg_D28_2_" --> "_4"
new_prefixes <- sapply(new_barcodes, GetPrefix) # 4 unique prefixes: "HuOrg_D16_1_" "HuOrg_D16_2_" "HuOrg_D28_1_" "HuOrg_D28_2_"
new_suffixes <- sapply(new_barcodes, GetSuffix) # 1 unique suffixes, ignore
new_cores <- sapply(new_barcodes, GetBarcodeCore)
old_prefixes <- sapply(old_barcodes, GetPrefix) # 4 unique prefixes: "1_" "2_" "3_" "4_"

modified_barcodes <- sapply(seq_along(old_barcodes), function(i) {
  prefix <- GetPrefix(old_barcodes[i])
  core <- GetBarcodeCore(old_barcodes[i])
  # suffix <- GetSuffix(old_barcodes[i])
  if (prefix == "1_") {
    modified_prefix <- "HuOrg_D16_1_"
  } else if (prefix == "2_") {
    modified_prefix <- "HuOrg_D16_2_"
  } else if (prefix == "3_") {
    modified_prefix <- "HuOrg_D28_1_"
  } else if (prefix == "4_") {
    modified_prefix <- "HuOrg_D28_2_"
  }

  modified_barcode <- paste0(modified_prefix, core, "-1")
  print(paste("Modified barcode:", modified_barcode))
  stopifnot(modified_barcode %in% new_barcodes)
  print("Checking if modified barcode is in new_barcodes")
  # print(modified_barcode %in% new_barcodes)
  return(modified_barcode)
})

stopifnot(!any(duplicated(modified_barcodes))) # No duplicate

# With all the above tests passed, we can probably assume that `H15w_Z1_` --> 1_ and `H15w_Z2_` --> 2.
# Now map. 
barcodes_mapping <- cbind(old_barcodes, modified_barcodes)

merged_meta <- merge(old_meta_file, barcodes_mapping, by.x=c(barcode_slot), by.y=c("old_barcodes"))
merged_meta_path <- ConnectPath(metadata_dir, paste0(study_id, "_merged_meta.tsv"))
write.table(merged_meta, merged_meta_path, row.names = FALSE, sep="\t")

# GSE125449
study_id <- "GSE125449"
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
  print("JnJ metadata exists")
  old_meta_file <- read.csv(meta_path, sep ="\t")
  barcode_slot <- ifelse('BBrowser_barcodes' %in% colnames(old_meta_file), 'BBrowser_barcodes', "Barcodes")
  old_barcodes <- as.character(old_meta_file[[barcode_slot]])
} else { # Export metadata manually and read it
  print("JnJ metadata not exists, take barcodes.tsv from old study zip")
  old_barcodes <- as.character(read.csv(ConnectPath(old_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
}

# production_metadata <- read.csv("/mnt2/vu/script/metadata/GSE124310/GSE124310_metadata_1632386464074.tsv", sep="\t")
new_barcodes <- as.character(read.csv(ConnectPath(new_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)

new_barcodes[1:5] #"HuOrg_D16_1_AAACCTGAGCGTAATA-1"
old_barcodes[1:5] #1_AAACCTGAGCGTAATA

# Need to parse the barcode part and check if unique 
truncated_new_barcodes <- stringr::str_extract(new_barcodes, "[ACGT]{16}")
length(unique(truncated_new_barcodes)) == length(truncated_new_barcodes)
stopifnot(length(unique(old_barcodes)) == length(old_barcodes))

# Try to change prefix  "Set1_" --> "_1" and  "Set2_" --> "_2"
new_prefixes <- sapply(new_barcodes, GetPrefix) # 
print("Unique prefixes in new barcodes: ")
unique(new_prefixes)
new_suffixes <- sapply(new_barcodes, GetSuffix) # 
print("Unique suffixes in new barcodes: ")
unique(new_suffixes)
new_cores <- sapply(new_barcodes, GetBarcodeCore)

old_prefixes <- sapply(old_barcodes, GetPrefix) # 
print("Unique prefixes in old barcodes: ")
unique(old_prefixes)
old_suffixes <- sapply(old_barcodes, GetSuffix) # 
print("Unique suffixes in old barcodes: ")
unique(old_suffixes)

stopifnot(length(unique(old_prefixes)) == length(unique(new_prefixes)))
stopifnot(length(unique(old_suffixes)) == length(unique(new_suffixes)))

print("Check that the count of suffixes in new and old barcodes are the same, if not, they should at least be similar")
all(plyr::count(new_suffixes) == plyr::count(old_suffixes)) # --> If TRUE, we can just focus on prefix

# A lot of suffixes

modified_barcodes <- sapply(seq_along(old_barcodes), function(i) {
  prefix <- GetPrefix(old_barcodes[i])
  core <- GetBarcodeCore(old_barcodes[i])
  suffix <- GetSuffix(old_barcodes[i])
  if (prefix == "1_") {
    modified_prefix <- "Set1_"
  } else if (prefix == "2_") {
    modified_prefix <- "Set2_"
  }

  modified_barcode <- paste0(modified_prefix, core, suffix)
  print(paste("Modified barcode:", modified_barcode))
  stopifnot(modified_barcode %in% new_barcodes)
  print("Checking if modified barcode is in new_barcodes")
  # print(modified_barcode %in% new_barcodes)
  return(modified_barcode)
})

stopifnot(!any(duplicated(modified_barcodes))) # No duplicate

# With all the above tests passed, we can probably assume that `H15w_Z1_` --> 1_ and `H15w_Z2_` --> 2.
# Now map. 
barcodes_mapping <- cbind(old_barcodes, modified_barcodes)

merged_meta <- merge(old_meta_file, barcodes_mapping, by.x=c(barcode_slot), by.y=c("old_barcodes"))
merged_meta_path <- ConnectPath(metadata_dir, paste0(study_id, "_merged_meta.tsv"))
write.table(merged_meta, merged_meta_path, row.names = FALSE, sep="\t")



#GSE125527
study_id <- "GSE125527"
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
  print("JnJ metadata exists")
  old_meta_file <- read.csv(meta_path, sep ="\t")
  barcode_slot <- ifelse('BBrowser_barcodes' %in% colnames(old_meta_file), 'BBrowser_barcodes', "Barcodes")
  old_barcodes <- as.character(old_meta_file[[barcode_slot]])
} else { # Export metadata manually and read it
  print("JnJ metadata not exists, take barcodes.tsv from old study zip")
  old_barcodes <- as.character(read.csv(ConnectPath(old_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
}

# production_metadata <- read.csv("/mnt2/vu/script/metadata/GSE124310/GSE124310_metadata_1632386464074.tsv", sep="\t")
new_barcodes <- as.character(read.csv(ConnectPath(new_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)

new_barcodes[1:5] #"HuOrg_D16_1_AAACCTGAGCGTAATA-1"
old_barcodes[1:5] #1_AAACCTGAGCGTAATA

# Need to parse the barcode part and check if unique 
truncated_new_barcodes <- stringr::str_extract(new_barcodes, "[ACGT]{16}")
length(unique(truncated_new_barcodes)) == length(truncated_new_barcodes)
stopifnot(length(unique(old_barcodes)) == length(old_barcodes))

# Try to change prefix  "Set1_" --> "_1" and  "Set2_" --> "_2"
new_prefixes <- sapply(new_barcodes, GetPrefix) # 
print("Unique prefixes in new barcodes: ")
unique(new_prefixes)
new_suffixes <- sapply(new_barcodes, GetSuffix) # 
print("Unique suffixes in new barcodes: ")
unique(new_suffixes)
new_cores <- sapply(new_barcodes, GetBarcodeCore)

old_prefixes <- sapply(old_barcodes, GetPrefix) # 
print("Unique prefixes in old barcodes: ")
unique(old_prefixes)
old_suffixes <- sapply(old_barcodes, GetSuffix) # 
print("Unique suffixes in old barcodes: ")
unique(old_suffixes)

stopifnot(length(unique(old_prefixes)) == length(unique(new_prefixes)))
stopifnot(length(unique(old_suffixes)) == length(unique(new_suffixes)))

print("Check that the count of suffixes in new and old barcodes are the same, if not, they should at least be similar")
all(plyr::count(new_suffixes) == plyr::count(old_suffixes)) # --> If TRUE, we can just focus on prefix
# --> For this study, prefixes is not necessary --> TODO: WHAT? Remove prefixes from new_barcodes
# A lot of suffixes

modified_barcodes <- sapply(seq_along(new_barcodes), function(i) {
  prefix <- GetPrefix(new_barcodes[i])
  core <- GetBarcodeCore(new_barcodes[i])
  suffix <- GetSuffix(new_barcodes[i])
  
  modified_barcode <- paste0(core, suffix)
  print(paste("Modified barcode:", modified_barcode))
  stopifnot(modified_barcode %in% old_barcodes)
  print("Checking if modified barcode is in new_barcodes")
  # print(modified_barcode %in% new_barcodes)
  return(modified_barcode)
})

stopifnot(!any(duplicated(modified_barcodes))) # No duplicate
stopifnot(all(modified_barcodes %in% old_barcodes))
# With all the above tests passed, we can probably assume that `H15w_Z1_` --> 1_ and `H15w_Z2_` --> 2.
# Now map. 
barcodes_mapping <- cbind(new_barcodes, modified_barcodes)

merged_meta <- merge(old_meta_file, barcodes_mapping, by.x=c(barcode_slot), by.y=c("modified_barcodes"))
merged_meta_path <- ConnectPath(metadata_dir, paste0(study_id, "_merged_meta.tsv"))
write.table(merged_meta, merged_meta_path, row.names = FALSE, sep="\t")



# GSE125881
study_id <- "GSE125881"
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
  print("JnJ metadata exists")
  old_meta_file <- read.csv(meta_path, sep ="\t")
  barcode_slot <- ifelse('BBrowser_barcodes' %in% colnames(old_meta_file), 'BBrowser_barcodes', "Barcodes")
  old_barcodes <- as.character(old_meta_file[[barcode_slot]])
} else { # Export metadata manually and read it
  print("JnJ metadata not exists, take barcodes.tsv from old study zip")
  old_barcodes <- as.character(read.csv(ConnectPath(old_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)
}

# production_metadata <- read.csv("/mnt2/vu/script/metadata/GSE124310/GSE124310_metadata_1632386464074.tsv", sep="\t")
new_barcodes <- as.character(read.csv(ConnectPath(new_study, 'main', 'barcodes.tsv'), sep='\t', header=FALSE)$V1)

new_barcodes[1:5] #"HuOrg_D16_1_AAACCTGAGCGTAATA-1"
old_barcodes[1:5] #1_AAACCTGAGCGTAATA

# Need to parse the barcode part and check if unique 
truncated_new_barcodes <- stringr::str_extract(new_barcodes, "[ACGT]{16}")
length(unique(truncated_new_barcodes)) == length(truncated_new_barcodes)
stopifnot(length(unique(old_barcodes)) == length(old_barcodes))

# Try to change prefix  "Set1_" --> "_1" and  "Set2_" --> "_2"
new_prefixes <- sapply(new_barcodes, GetPrefix) # 
print("Unique prefixes in new barcodes: ")
unique(new_prefixes)
new_suffixes <- sapply(new_barcodes, GetSuffix) # 
print("Unique suffixes in new barcodes: ")
unique(new_suffixes)
new_cores <- sapply(new_barcodes, GetBarcodeCore)

old_prefixes <- sapply(old_barcodes, GetPrefix) # 
print("Unique prefixes in old barcodes: ")
unique(old_prefixes)
old_suffixes <- sapply(old_barcodes, GetSuffix) # 
print("Unique suffixes in old barcodes: ")
unique(old_suffixes)

## Seems like the suffix of old barcodes is the same as prefix of new_barcodes
print("Check that the count of suffixes in new and old barcodes are the same, if not, they should at least be similar")

new_prefixes_count <- plyr::count(new_prefixes)
colnames(new_prefixes_count)[1] <- 'new_prefix'

old_suffixes_count <- plyr::count(old_suffixes)
colnames(old_suffixes_count)[1] <- 'old_suffix'

counts <- cbind(old_suffixes_count, new_prefixes_count)
stopifnot(plyr::count(new_prefixes)$freq == plyr::count(old_suffixes)$freq)

unique_old_suffixes <- as.character(counts$old_suffix)
unique_new_prefixes <- as.character(counts$new_prefix)
# --> Rau ong no chap duoi pa` kia
modified_barcodes <- sapply(seq_along(new_barcodes), function(i) {
  prefix <- GetPrefix(new_barcodes[i])
  core <- GetBarcodeCore(new_barcodes[i])
  j <- match(prefix, unique_new_prefixes)
  suffix <- unique_old_suffixes[j]

  # Checkpoint
  trimmed_prefixed <- stringi::stri_trim_both(prefix, pattern = "[-_]", negate=TRUE)    
  trimmed_suffixed <- stringi::stri_trim_both(suffix, pattern = "[-_]", negate=TRUE)  
  stopifnot(trimmed_prefixed == trimmed_suffixed)
  #

  modified_barcode <- paste0(core, suffix)
  print(paste("Modified barcode:", modified_barcode))
  stopifnot(modified_barcode %in% old_barcodes)
  print("Checking if modified barcode is in new_barcodes")
  # print(modified_barcode %in% new_barcodes)
  return(modified_barcode)
})

stopifnot(!any(duplicated(modified_barcodes))) # No duplicate
stopifnot(all(modified_barcodes %in% old_barcodes))

barcodes_mapping <- cbind(new_barcodes, modified_barcodes)
print("Sample barcode_mapping...")
print(barcodes_mapping[10000:10005,])

merged_meta <- merge(old_meta_file, barcodes_mapping, by.x=c(barcode_slot), by.y=c("modified_barcodes"))
merged_meta_path <- ConnectPath(metadata_dir, paste0(study_id, "_merged_meta.tsv"))
write.table(merged_meta, merged_meta_path, row.names = FALSE, sep="\t")


# GSE125970
study_id <- "GSE125970"
OUTPUT_DIR <- "/mnt2/vu/script/output"
META_DIR <- "/mnt1/tracy/jnj/metadata"

res <- CheckStudyMetadata(study_id)
old_barcodes <- res$old_barcodes
new_barcodes <- res$new_barcodes
# Need to parse the barcode part and check if unique 
CheckCoreUniqueness(new_barcodes) # Should be TRUE

diagnosis <- DiagnoseMisMatch(new_barcodes, old_barcodes)
diagnosis$res[[diagnosis$new]][1:5]
diagnosis$res[[diagnosis$old]][1:5]
counts <- diagnosis$counts
# Seems like "Ileum-1_"  "Ileum-2_"  "Colon-1_"  "Colon-2_"  "Rectum-1_" "Rectum-2_" 
# correspond to "1_" "2_" "3_" "4_" "5_" "6_", Not sure about exact correspondence though --> Check with count

# Try to change prefix
#  "Colon-1_" --> "6_" 
#  "Colon-2_" --> "_5"
#  "Ileum-1_" --> "4_"
#  "Ileum-2_" --> "3_"
#  "Rectum-1_" --> "2_"
#  "Rectum-2_" --> "1_"

new_terms <- as.character(counts$x1)
old_terms <- as.character(counts$x2)

new_cores <- sapply(new_barcodes, GetBarcodeCore)
old_cores <- sapply(old_barcodes, GetBarcodeCore)
modified_barcodes <- sapply(seq_along(old_barcodes), function(i) {
  part <- ParseBarcode(old_barcodes[i], diagnosis$old)
  core <- GetBarcodeCore(old_barcodes[i])

  j <- match(prefix, unique_new_prefixes)
  suffix <- unique_old_suffixes[j]

  # Checkpoint
  trimmed_prefixed <- stringi::stri_trim_both(prefix, pattern = "[-_]", negate=TRUE)    
  trimmed_suffixed <- stringi::stri_trim_both(suffix, pattern = "[-_]", negate=TRUE)  
  stopifnot(trimmed_prefixed == trimmed_suffixed)
  #

  modified_barcode <- paste0(core, suffix)
  print(paste("Modified barcode:", modified_barcode))
  stopifnot(modified_barcode %in% old_barcodes)
  print("Checking if modified barcode is in new_barcodes")
  # print(modified_barcode %in% new_barcodes)
  return(modified_barcode)
})

stopifnot(!any(duplicated(modified_barcodes))) # No duplicate
stopifnot(all(modified_barcodes %in% old_barcodes))

barcodes_mapping <- cbind(new_barcodes, modified_barcodes)
print("Sample barcode_mapping...")
print(barcodes_mapping[10000:10005,])

merged_meta <- merge(old_meta_file, barcodes_mapping, by.x=c(barcode_slot), by.y=c("modified_barcodes"))
merged_meta_path <- ConnectPath(metadata_dir, paste0(study_id, "_merged_meta.tsv"))
write.table(merged_meta, merged_meta_path, row.names = FALSE, sep="\t")









### REMAP WITH METADATA FROM TRANG'S ROUND 1
batch_4 = "/mnt2/bioturing_data/files/batch_4.tsv"
study_list = as.character(read.table(file = batch_4, sep = '\t', header = FALSE)$V1) 
meta_dir = "/mnt1/tracy/index/meta_map/output"
output_dir <- "/mnt2/vu/script/output"
devtools::load_all('../NoraSC')

for (study_id in study_list) {
  print(study_id)
  meta_path <- ConnectPath(meta_dir, paste0(study_id, ".tsv"))
  stopifnot(file.exists(meta_path))
  metadata <- read.csv(meta_path, sep="\t")
  stopifnot("Barcodes_1st" %in% colnames(metadata))
  print(colnames(metadata)[1])
}

# Every study in batch_4 has "Barcodes_1st", now compare that columns with hdf5's barcodes
for (study_id in study_list) {
  print(study_id)
  res <- CheckStudyMetadata(study_id)
  old_barcodes <- res$old_barcodes
  new_barcodes <- res$new_barcodes
  
  print(paste("Number of old and new barcodes:", length(old_barcodes), length(new_barcodes)))
  new_cores <- sapply(new_barcodes, GetBarcodeCore)
  old_cores <- sapply(old_barcodes, GetBarcodeCore)

  print("CHECK if All old_barcodes core are in new barcodes")
  print(all(old_cores %in% new_cores))
}

