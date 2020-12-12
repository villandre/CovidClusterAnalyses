# There's a long site of potential homoplasic sites provided with Van Dorp 2020.

# homoplasicSitesVanDorp <- read.csv("data/VanDorpHomoplasies.csv")$bp
# intersect(c(mainHomoplasicSites, otherHomoplasicSites), homoplasicSitesVanDorp) # 29 out of 30 sites proposed as homoplasic by De Maio were also proposed as homoplasic by Van Dorp.
# Function returns the name of the file where the aligned sequences were saved.

alignAndSaveCanadianData <- function(SARScov2reference, folderForSequences, folderWhereToSaveTheResult, patternForImport = ".pass.fasta", patternForSeqName = "(?<=/)(L|MCG|MGC|CHAL-|HSU-|JUS-|CHUM-|HCLM-|HDS-|HGA-|HMR-|HSE-|HVE-|JEW-|RIM-|S).+?(?=\\.)", numThreads, extremitiesBoundaries = c(55, 29804), homoplasicSites = c(c(187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408, 14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681, 28077, 28826, 28854, 29700), c(4050, 13402, 11083, 15324, 21575)), clustalExecutableName = "clustalo-1.2.4-Ubuntu-x86_64", resolutionRequirement = 0.95, numSites = 29903) {
  refLength <- length(SARScov2reference[[1]])

  filesToImport <- list.files(path = folderForSequences, pattern = patternForImport, recursive = TRUE, full.names = TRUE)

  # Do we have repeated sequences?
  seqNames <- stringr::str_extract(filesToImport, pattern = patternForSeqName)
  nonStandardSeqsPos <- which(is.na(seqNames))
  seqNames <- seqNames[-nonStandardSeqsPos]

  # We prioritise nanopore -> illumina -> mgi
  uniqueFilenames <- sapply(split(filesToImport[-nonStandardSeqsPos], seqNames), function(x) {
    returnValue <- x
    if (length(x) > 1) {
      returnValue <- stringr::str_subset(x, pattern = "nanopore")
      if (length(returnValue) == 0) returnValue <- stringr::str_subset(x, pattern = "illumina")
    }
    returnValue
  })

  sequencesList <- lapply(uniqueFilenames, ape::read.FASTA)

  resolutionThreshold <- resolutionRequirement * refLength

  includeIndex <- sapply(sequencesList, function(sequenceDNAbin) {
    if (length(sequenceDNAbin[[1]]) < resolutionThreshold) return(FALSE)

    nucleoTable <- table(as.character(sequenceDNAbin[1]))
    if (!("n" %in% names(nucleoTable))) return(TRUE)
    if (nucleoTable[["n"]] < refLength - resolutionThreshold) return(TRUE)
    return(FALSE)
  })

  funForLapply <- function(seqIndex, sequencesList, numSites, clustalExecutableName, SARScov2reference) {
    cat("Processing sequence ", seqIndex, "... ")
    alignedSeq <- sequencesList[[seqIndex]]
    if (length(sequencesList[[seqIndex]][[1]]) != numSites) {
      alignedSeq <- tryCatch(expr = ape::clustal(x = sequencesList[[seqIndex]], y = SARScov2reference,  exec = clustalExecutableName), error = function(e) e)
      if ("DNAbin" %in% class(alignedSeq)) {
        alignedSeq <- alignedSeq[1, which(as.character(alignedSeq[2, ]) != "-")]
      }
    }
    cat("Done! \n\n")
    alignedSeq
  }
  sequences <- NULL
  if (numThreads > 1) {
    cl <- parallel::makeForkCluster(numThreads)
    sequences <- parallel::parLapply(which(includeIndex), cl = cl, fun = funForLapply, sequencesList = sequencesList, numSites = numSites, clustalExecutableName = clustalExecutableName, SARScov2reference = SARScov2reference)
    parallel::stopCluster(cl)
  } else {
    sequences <- lapply(which(includeIndex), funForLapply, sequencesList = sequencesList, numSites = numSites, clustalExecutableName = clustalExecutableName, SARScov2reference = SARScov2reference)
  }
  SARScov2dataAlignedCanada <- do.call(rbind, lapply(sequences, as.matrix))

  # We got the following list of sites to exclude from De Maio 2020. They are considered sites where data is of lower quality (sequence extremities) or sites with homoplasies.

  extremitiesToExclude <- c(1:extremitiesBoundaries[[1]], extremitiesBoundaries[[2]]:29903)
  sitesToKeep <- setdiff(1:refLength, c(extremitiesToExclude, homoplasicSites))
  currentDateString <- as.character(Sys.Date())
  currentDateStringCorrected <- stringr::str_replace_all(currentDateString, pattern = "-", replacement = "")
  sequencesFilename <- paste(folderWhereToSaveTheResult, "/SARScov2dataAlignedCanada_", currentDateStringCorrected, ".fasta", sep = "")
  ape::write.dna(SARScov2dataAlignedCanada[ , sitesToKeep] , file = sequencesFilename, format = "fasta", nbcol = 1, colw = 100)
  cat("Saved sequences in ", sequencesFilename, "\n")
  sequencesFilename
}

extractMetadata <- function(DNAbinObject, folderForMetadata, patternForMetadataFiles = "minimal", patternInSequenceNames = "(?<=(c|C)anada/Qc-).+(?=/2020)", seqNameColumn = "sample", sampleDateColumnName = "sample_date") {
  DNAbinObject <- as.matrix(DNAbinObject)
  sequencesNames <- stringr::str_extract(rownames(DNAbinObject), patternInSequenceNames)
  rownames(DNAbinObject) <- sequencesNames

  metadataFiles <- list.files(path = folderForMetadata, pattern = patternForMetadataFiles, full.names = TRUE)
  metadataList <- lapply(metadataFiles, read.table, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  metadata <- unique(do.call(rbind, metadataList))

  # Not all sequences have a precise date. We'll only consider sequences with a precise sampling date.

  metadataWithCompleteDate <- metadata[nchar(metadata[ , sampleDateColumnName]) > 7, ]
  metadataWithCompleteDate[ , sampleDateColumnName] <- as.POSIXct(metadataWithCompleteDate[ , sampleDateColumnName])

  funForLapply <- function(seqName) {
    rowNumber <- stringr::str_which(string = metadataWithCompleteDate[ , seqNameColumn], pattern = seqName)
    if (length(rowNumber) > 1) {
      rowNumber <- rowNumber[which.min(metadataWithCompleteDate$date[rowNumber])]
    }
    returnValue <- metadataWithCompleteDate[rowNumber, ]
    if (length(rowNumber) == 0) {
      returnValue <- c(strain = seqName, date = NA, country = "Canada", division = "Quebec", rss = NA)
    }
    returnValue
  }
  matchingMetadataRowsList <- lapply(sequencesNames, funForLapply)
  matchingMetadataRows <- do.call(rbind, matchingMetadataRowsList)

  sequencesToKeepIndices <- which(!is.na(matchingMetadataRows$date))
  sequencesDates <- matchingMetadataRows$date[sequencesToKeepIndices]
  sequencesRegions <- matchingMetadataRows$division[sequencesToKeepIndices]
  names(sequencesDates) <- names(sequencesRegions) <- sequencesNames[sequencesToKeepIndices]
  list(timestamps = sequencesDates, regionStamps = sequencesRegions)
}
