# There's a long site of potential homoplasic sites provided with Van Dorp 2020.

# homoplasicSitesVanDorp <- read.csv("data/VanDorpHomoplasies.csv")$bp
# intersect(c(mainHomoplasicSites, otherHomoplasicSites), homoplasicSitesVanDorp) # 29 out of 30 sites proposed as homoplasic by De Maio were also proposed as homoplasic by Van Dorp.
# Function returns the name of the file where the aligned sequences were saved.

alignAndSaveCanadianData <- function(SARScov2reference, folderForSequences, folderWhereToSaveTheResult, patternForImport = ".pass.fasta", patternForSeqName = "(?<=/)(L|MCG|MGC|CHAL-|HSU-|JUS-|CHUM-|HCLM-|HDS-|HGA-|HMR-|HSE-|HVE-|JEW-|RIM-|S).+?(?=\\.)", numThreads, extremitiesBoundaries = c(55, 29804), homoplasicSites = c(c(187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408, 14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681, 28077, 28826, 28854, 29700), c(4050, 13402, 11083, 15324, 21575)), resolutionRequirement = 0.95, numSites = 29903) {
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

  funForLapply <- function(seqIndex, sequencesList, numSites, SARScov2reference) {
    cat("Processing sequence ", seqIndex, "... ")
    alignedSeq <- sequencesList[[seqIndex]]
    # if (length(sequencesList[[seqIndex]][[1]]) != numSites) { # No chance to take if there's a sequence with one missing nucleotide on one side, and one added on the other side.
    alignedSeq <- tryCatch(expr = ape::clustalomega(x = c(sequencesList[[seqIndex]], SARScov2reference)), error = function(e) e)
    if ("DNAbin" %in% class(alignedSeq)) {
      alignedSeq <- alignedSeq[1, which(as.character(alignedSeq[2, ]) != "-")]
    }
    # }
    cat("Done! \n\n")
    alignedSeq
  }
  sequences <- NULL
  if (numThreads > 1) {
    cl <- parallel::makeForkCluster(numThreads)
    sequences <- parallel::parLapply(which(includeIndex), cl = cl, fun = funForLapply, sequencesList = sequencesList, numSites = numSites, SARScov2reference = SARScov2reference)
    parallel::stopCluster(cl)
  } else {
    sequences <- lapply(which(includeIndex), funForLapply, sequencesList = sequencesList, numSites = numSites, SARScov2reference = SARScov2reference)
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

extractAndSaveCanadianMetadata <- function(DNAbinObject, folderForMetadata, folderToSaveResult, patternForMetadataFiles = "minimal", patternInSequenceNames = "(?<=(c|C)anada/Qc-).+(?=/2020)", seqNameColumn = "sample", sampleDateColumnName = "sample_date", regionColumnName = "division") {
  DNAbinObject <- as.matrix(DNAbinObject)
  # sequencesNames <- stringr::str_extract(rownames(DNAbinObject), patternInSequenceNames)
  # rownames(DNAbinObject) <- sequencesNames

  metadataFiles <- list.files(path = folderForMetadata, pattern = patternForMetadataFiles, full.names = TRUE)
  metadataList <- lapply(metadataFiles, read.table, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  colnamesToKeep <- Reduce(f = "intersect", sapply(metadataList, colnames))
  metadata <- unique(do.call(rbind, lapply(metadataList, function(listElement) listElement[ , colnamesToKeep])))

  # Not all sequences have a precise date. We'll only consider sequences with a precise sampling date.

  metadataWithCompleteDate <- metadata[nchar(metadata[ , sampleDateColumnName]) > 7, ]
  metadataWithCompleteDate[ , sampleDateColumnName] <- as.POSIXct(metadataWithCompleteDate[ , sampleDateColumnName])
  metadataIndicesToKeep <- which(!duplicated(metadataWithCompleteDate[ , seqNameColumn]))

  seqIndicesMatchingNamesInMetadata <- lapply(metadataWithCompleteDate[metadataIndicesToKeep, seqNameColumn], FUN = function(seqName) {
    matchNumbers <- grep(pattern = seqName, x = rownames(DNAbinObject))
    if (length(matchNumbers) == 0) matchNumbers <- NA
    matchNumbers[[1]] # There might be repeated rows...
  }) # Some sequences mentioned in metadata have been flagged or rejected. What does it mean?
  seqIndicesToKeepWithNAs <- do.call("c", seqIndicesMatchingNamesInMetadata)
  sequencesToKeepIndices <- seqIndicesToKeepWithNAs[!is.na(seqIndicesToKeepWithNAs)]
  sequencesDates <- metadataWithCompleteDate[metadataIndicesToKeep[!is.na(seqIndicesToKeepWithNAs)], sampleDateColumnName]
  sequencesRegions <- metadataWithCompleteDate[metadataIndicesToKeep[!is.na(seqIndicesToKeepWithNAs)], regionColumnName]
  names(sequencesDates) <- names(sequencesRegions) <- metadataWithCompleteDate[metadataIndicesToKeep[!is.na(seqIndicesToKeepWithNAs)], seqNameColumn]
  # Some timestamps are wrong. We must remove sequences with faulty timestamps.
  wrongTimestamps <- sequencesDates < as.POSIXct("2020-01-01 UTC")
  sequencesToKeepIndices <- sequencesToKeepIndices[!wrongTimestamps]
  sequencesRegions <- sequencesRegions[!wrongTimestamps]
  sequencesDates <- sequencesDates[!wrongTimestamps]
  listToSave <- list(timestamps = sequencesDates[order(sequencesToKeepIndices)], regionStamps = tolower(sequencesRegions)[order(sequencesToKeepIndices)], seqsToKeepIndices = sort(sequencesToKeepIndices))
  currentDateString <- as.character(Sys.Date())
  currentDateStringCorrected <- stringr::str_replace_all(currentDateString, pattern = "-", replacement = "")
  filename <- paste(folderToSaveResult, "/canadianMetadata_", currentDateStringCorrected, ".rds", sep = "")
  saveRDS(listToSave, file = filename, compress = TRUE)
  cat("Canadian metadata saved in", filename, "\n")
  filename
}

alignAndSaveGISAIDdata <- function(SARScov2reference, folderForSequences, folderWhereToSaveTheResult, GISAIDfastaFilename, numThreads, extremitiesBoundaries = c(55, 29804), homoplasicSites = c(c(187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408, 14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681, 28077, 28826, 28854, 29700), c(4050, 13402, 11083, 15324, 21575)), resolutionRequirement = 0.95, numSites = 29903, aligned = FALSE) {
  GISAIDdata <- ape::read.FASTA(GISAIDfastaFilename)

  refLength <- length(SARScov2reference[[1]])

  resolutionThreshold <- resolutionRequirement * refLength

  includeIndex <- sapply(seq_along(GISAIDdata), function(index) {
    if (length(GISAIDdata[[index]]) < resolutionThreshold) return(FALSE)

    nucleoTable <- table(as.character(GISAIDdata[index]))
    if (!("n" %in% names(nucleoTable))) return(TRUE)
    if (nucleoTable[["n"]] < refLength - resolutionThreshold) return(TRUE)
    return(FALSE)
  })

  seqsToAlignIndices <- which(includeIndex)

  funForLapply <- function(SARScov2seqIndex, SARScov2reference, SARScov2dataV0) {
    alignedPair <- ape::clustalomega(x = c(SARScov2dataV0[SARScov2seqIndex], SARScov2reference))
    # alignedSeqs[1, which(as.character(alignedSeqs[2, ]) != "-")]
    alignedPair[1, ]
  }

  alignedSeqs <- NULL
  if (!aligned) {
  if (numThreads > 1) {
    cl <- parallel::makeForkCluster(numThreads)
    alignedSeqs <- parallel::parLapply(seqsToAlignIndices, cl = cl, fun = funForLapply, SARScov2reference = SARScov2reference, SARScov2dataV0 = GISAIDdata) ## Not all sequences have the same length! How do we deal with potential insertions? In practice, what's the difference between "n" and "-"? There might be an insertion in hCoV-19/USA/WA13-UW9/2020|EPI_ISL_413601|2020-03-02 (Seq. 92).
  parallel::stopCluster(cl)
  } else {
    alignedSeqs <- lapply(seqsToAlignIndices, funForLapply, SARScov2reference = SARScov2reference, SARScov2dataV0 = GISAIDdata)
  }
  GISAIDdata <- Reduce(f = "c", alignedSeqs)
  }

  extremitiesToExclude <- c(1:55, 29804:29903)
  mainHomoplasicSites <- c(187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408, 14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681, 28077, 28826, 28854, 29700)
  otherHomoplasicSites <- c(4050, 13402, 11083, 15324, 21575)

  sitesToKeep <- setdiff(1:refLength, c(extremitiesToExclude, mainHomoplasicSites, otherHomoplasicSites))

  GISAIDdata <- as.matrix(GISAIDdata)[ , sitesToKeep]

  currentDateString <- as.character(Sys.Date())
  currentDateStringCorrected <- stringr::str_replace_all(currentDateString, pattern = "-", replacement = "")
  sequencesFilename <- paste(folderWhereToSaveTheResult, "/SARScov2dataAlignedGISAIDclean_", currentDateStringCorrected, ".fasta", sep = "")

  ape::write.dna(GISAIDdata, file = sequencesFilename, format = "fasta", nbcol = 1, colw = 100)
  cat("Saved GISAID sequences in", sequencesFilename, "\n")
  sequencesFilename
}

extractAndSaveGISAIDmetadata <- function(GISAIDdata, regionsToSampleFrom = c("canada", "australia", "china", "belgium", "denmark", "france",  "england", "finland", "germany", "italy", "portugal", "russia", "spain", "sweden", "turkey", "usa", "wales"), seed, numToSamplePerRegion = 10, folderToSaveResult, fixRegionLabelsList = list(list(replacement = "china", original = c("anhui", "beijing", "chongqing", "fujian", "fuyang", "ganzhou", "guangdong", "guangzhou", "hangzhou", "hefei", "henan", "hong kong", "jian", "jiangsu", "jiangxi", "jingzhou", "jiujiang", "lishui", "nanchang", "pingxiang", "shandong", "shanghai", "shangrao", "shenzhen", "sichuan", "tianmen", "wuhan", "yunnan", "zhejiang")))) {
  rawSeqNames <- rownames(GISAIDdata)
  if (is.null(rawSeqNames)) {
    rawSeqNames <- names(GISAIDdata)
  }
  samplingDates <- stringr::str_extract(rawSeqNames, "(?<=\\|)2020.*$")
  regionStampsFull <- stringr::str_to_lower(stringr::str_extract(rawSeqNames, pattern = "(?<=hCoV-19/)[a-zA-Z\\s]+(?=/)"))
  for (i in 1:length(fixRegionLabelsList)) {
    regionStampsFull <- replace(regionStampsFull, list = which(regionStampsFull %in% fixRegionLabelsList[[i]]$original), values = fixRegionLabelsList[[i]]$replacement)
  }

  seqsToSampleFromRegionCrit <- regionStampsFull %in% regionsToSampleFrom
  seqsToSampleFromTimeCrit <- nchar(samplingDates) == 10
  seqsToSampleFrom <- which(seqsToSampleFromRegionCrit & seqsToSampleFromTimeCrit)

  set.seed(seed)
  seqIndicesByRegion <- tapply(seqsToSampleFrom, INDEX = regionStampsFull[seqsToSampleFrom], FUN = function(seqIndices) {
    numToSample <- min(numToSamplePerRegion, length(seqIndices))
    sample(seqIndices, size = numToSample, replace = FALSE)
  }, simplify = FALSE)
  seqsToKeepIndices <- sort(do.call("c", seqIndicesByRegion))
  timestampsPOSIXct <- as.POSIXct(samplingDates[seqsToKeepIndices]) # Time in days
  regionStamps <- regionStampsFull[seqsToKeepIndices]
  names(timestampsPOSIXct) <- names(regionStamps) <- rawSeqNames[seqsToKeepIndices]
  listToSave <- list(timestamps = timestampsPOSIXct, regionStamps = tolower(regionStamps), seqsToKeepIndices = seqsToKeepIndices)
  currentDateString <- as.character(Sys.Date())
  currentDateStringCorrected <- stringr::str_replace_all(currentDateString, pattern = "-", replacement = "")
  filename <- paste(folderToSaveResult, "/GISAIDmetadata_", currentDateStringCorrected, ".rds", sep = "")
  saveRDS(listToSave, file = filename, compress = TRUE)
  cat("GISAID metadata saved in", filename, "\n")
  filename
}
