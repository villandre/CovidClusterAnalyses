plotClustersTile <- function(covidClusterObject, minClusSize = 10, textScaleFactor = 1, oneFilePerCluster = FALSE, outputFolder, device = "pdf", numClusters = Inf, clustersToPlot = "hier", controlForDevice = NULL) {
  plotObject <- covidClusterObject$objectToPlot
  clusterIndices <- covidClusterObject$MAPclusters
  if (clustersToPlot == "hier") {
    clusterIndices <- covidClusterObject$hierarchicalClusters
  }
  clusSizes <- table(clusterIndices)
  largeClusters <- as.numeric(names(clusSizes))[clusSizes >= minClusSize]
  largeClusters <- head(largeClusters, n = numClusters)

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
