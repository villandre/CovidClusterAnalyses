plotClustersTile <- function(covidClusterObject, minClusSize = 10, textScaleFactor = 1, oneFilePerCluster = FALSE, outputFolder, device = "pdf", numClusters = Inf, clustersToPlot = "hier", controlForDevice = NULL) {
  seqNames <- rownames(covidClusterObject$adjMatrix)
  frameToPlot <- as.data.frame(Matrix::mat2triplet(covidClusterObject$adjMatrix))
  colnames(frameToPlot) <- c("Sequence_x", "Sequence_y", "CoclusteringRate")
  frameToPlot$Sequence_x <- factor(frameToPlot$Sequence_x, levels = seq_along(seqNames), labels = seqNames)
  frameToPlot$Sequence_y <- factor(frameToPlot$Sequence_y, levels = seq_along(seqNames), labels = seqNames)
  frameToPlot <- subset(frameToPlot, subset = CoclusteringRate > 1e-8)
  freqTable <- table(frameToPlot$Sequence_x)
  singletons <- names(freqTable)[which(freqTable == 1)]
  frameToPlot <- subset(frameToPlot, subset = !(Sequence_x %in% singletons))

  clusterIndices <- covidClusterObject$MAPclusters
  if (clustersToPlot == "hier") {
    clusterIndices <- covidClusterObject$hierarchicalClusters
  }
  clusSizes <- table(clusterIndices)
  largeClusters <- as.numeric(names(clusSizes))[clusSizes >= minClusSize]
  largeClusters <- head(largeClusters, n = numClusters)
  plotObject <-
    ggplot2::ggplot(data = frameToPlot, ggplot2::aes(x = Sequence_x, y = Sequence_y, fill = CoclusteringRate)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_steps(low = "white", high = "red", limit = c(0,1.01), space = "Lab", name = "Coclustering\nrate") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
    ggplot2::coord_fixed()

  if (oneFilePerCluster) {
    tilePlotList <- lapply(largeClusters, FUN = function(clusIndex) {
      seqsInClus <- names(clusterIndices)[clusterIndices == clusIndex]
      newPlotObject <- plotObject + ggplot2::scale_x_discrete(limits = seqsInClus) + ggplot2::scale_y_discrete(limits = seqsInClus)
      newPlotObject
    })
    filename <- paste(outputFolder, "/clustersTilePlots.", device, sep = "")
    argsForDevice <- c(list(file = filename), controlForDevice)
    do.call(device, args = argsForDevice) # Summons the plotting device
    do.call(gridExtra::grid.arrange, tilePlotList)
    dev.off()
    returnObject <- tilePlotList
  } else {
    seqsToKeep <- names(clusterIndices)[clusterIndices %in% largeClusters]
    plotObject <- plotObject + ggplot2::scale_x_discrete(limits = seqsToKeep) +  ggplot2::scale_y_discrete(limits = seqsToKeep) + ggplot2::theme(axis.text = ggplot2::element_text(colour = "blue", size = textScaleFactor))
    returnObject <- plotObject
    filename <- paste(outputFolder, "/clustersTilePlot.", device, sep = "")
    do.call(ggplot2::ggsave, args = list(filename = filename, plot = plotObject, device = device, controlForDevice))
  }
  returnObject
}

sparseMatrixFindBlockStartRows <- function(sparseMat) {
  triplets <- Matrix::mat2triplet(sparseMat)
  nonZeroIndices <- cbind(triplets$i, triplets$j)
  currentPos <- 1
  outputVector <- 1
  while (currentPos <= nrow(sparseMat)) {
    matrixColumn <- nonZeroIndices[nonZeroIndices$j == currentPos]
    diffsInCol <- diff(sort(matrixColumn))
    partOfBlock <- diffsInCol == 1
    currentPos <- match(FALSE, partOfBlock) + 2
    outputVector <- c(outputVector, currentPos) # Not memory-efficient, but probably not time-consuming
    if (is.na(currentPos)) currentPos <- length(diffsInCol) + 2
  }
  head(outputVector, n = -1) # The last element of outputVec is nrow(sparseMat) + 1, hence the call to head.
}
