MatchBarcodes <- function(bc.1, bc.2) {
  TrimMultiplexingTags <- function(bc.raw) {
    bc <- sapply(strsplit(bc.raw, "-"), function(x) {
      return(paste(x[1: min(1, length(x) - 1)], collapse = "-"))
    })
    return(bc)
  }

  TrimBioTuringPrefixes <- function(bc.raw) {
    bc <- sapply(strsplit(bc.raw, "_"), function(x) {
      return(paste(x[min(length(x), 2): length(x)], collapse = "_"))
    })
    return(bc)
  }

  type.1 <- GetBarcodesType(bc.1)
  type.2 <- GetBarcodesType(bc.2)

  # Handle multiplexing tags
  if (type.1$isMultiplexed != type.2$isMultiplexed) {
    if (type.1$isMultiplexed) {
      bc.1 <- TrimMultiplexingTags(bc.1)
    } else {
      bc.2 <- TrimMultiplexingTags(bc.2)
    }
  }

  # Handle BioTuring prefixes
  if (type.1$isBioTuring != type.2$isBioTuring) {
    if (type.1$isBioTuring) {
      bc.1 <- TrimBioTuringPrefixes(bc.1)
    } else {
      bc.2 <- TrimBioTuringPrefixes(bc.2)
    }
  }

  return(match(bc.1, bc.2))
}

GetBarcodesType <- function(bc.raw) {
  bc <- bc.raw[1: min(length(bc.raw), 1000)]
  type <- list(
    isBioTuring = all(grepl("^\\d+\\_.+", bc)),
    isMultiplexed = all(grepl("\\-\\d$", bc))
  )
  return(type)
}
