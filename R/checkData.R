identifyNonOutlierInternationalSeqs <- function(SARScov2reference, canadianDNAbinData, internationalDNAbinData, canadianIndices, internationalIndices, extremitiesBoundaries = c(55, 29804), homoplasicSites = c(c(187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408, 14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681, 28077, 28826, 28854, 29700), c(4050, 13402, 11083, 15324, 21575)), excessDistThreshold = 2) {
  SARScov2reference <- as.matrix(SARScov2reference)[ , -c(1:extremitiesBoundaries[[1]], extremitiesBoundaries[[2]]:29903, homoplasicSites)]
  distToRef <- lapply(
    list(Canada = as.matrix(canadianDNAbinData)[canadianIndices, ], International = as.matrix(internationalDNAbinData)[internationalIndices, ]), function(dataset) {
      distToRef <- sapply(1:nrow(dataset), function(seqIndex) {
        seqPair <- rbind(dataset[seqIndex, ], SARScov2reference[1, ])
        ape::dist.dna(seqPair, pairwise.deletion = TRUE, as.matrix = TRUE)[1, 2]
      })
    })
  print(lapply(distToRef, summary))
  cat("Maximum distance for Canadian sequences:", max(distToRef$Canada), "reached for sequence", names(canadianDNAbinData)[which.max(distToRef$Canada)], "\n")
  cat("Maximum distance for international sequences:", max(distToRef$International), "reached for sequence", names(internationalDNAbinData)[which.max(distToRef$International)], "\n")
  maxAllowedDistance <- max(distToRef$Canada) * excessDistThreshold
  internationalIndices[distToRef$International <= maxAllowedDistance]
} # internationalSeqsToKeepIndices <- identifyNonOutlierInternationalSeqs(SARScov2reference, alignedCanadianSequences, alignedGISAIDdata, canadianIndices = canadianMetadata$seqsToKeepIndices, internationalIndices = internationalMetadata$seqsToKeepIndices, excessDistThreshold = 2)

