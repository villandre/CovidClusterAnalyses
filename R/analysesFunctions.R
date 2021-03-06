plotClustersTile <- function(covidClusterObject, clusterIndices, minClusSize = 10, textScaleFactor = 1, oneFilePerCluster = FALSE, outputFolder, device = "pdf", numClusters = Inf, controlListGgsave = NULL) {
  seqNames <- names(covidClusterObject$MAPclusters)
  frameToPlot <- as.data.frame(Matrix::mat2triplet(Matrix::tril(covidClusterObject$adjMatrix, k = -1)))
  colnames(frameToPlot) <- c("Sequence_x", "Sequence_y", "CoclusteringRate")
  matchForX <- nrow(covidClusterObject$adjMatrix):1
  frameToPlot$Sequence_x <- matchForX[frameToPlot$Sequence_x]
  frameToPlot$Sequence_x <- factor(frameToPlot$Sequence_x, levels = seq_along(seqNames), labels = rev(seqNames))
  frameToPlot$Sequence_y <- factor(frameToPlot$Sequence_y, levels = seq_along(seqNames), labels = seqNames)

  clusSizes <- table(clusterIndices)
  largeClusters <- as.numeric(names(clusSizes))[clusSizes >= minClusSize]
  largeClusters <- head(largeClusters, n = numClusters)

  if (oneFilePerCluster) {
    produceOneClusTilePlot <- function(clusIndex) {
      seqsInClus <- names(clusterIndices)[clusterIndices == clusIndex]
      frameToPlot <- subset(frameToPlot, (as.character(Sequence_x) %in% seqsInClus) & (as.character(Sequence_y) %in% seqsInClus))
      frameToPlot$Sequence_x <- droplevels(frameToPlot$Sequence_x)
      frameToPlot$Sequence_y <- droplevels(frameToPlot$Sequence_y)
      # Defining plotObject outside the function and subsetting the rows/columns of the tile plot with scale_*_discrete(limits = ...) results in a glitch I still can't explain, where the lowermost tile ends up above the diagonal instead of below, like the rest.
      plotObject <-
        ggplot2::ggplot(data = frameToPlot, ggplot2::aes(x = Sequence_y, y = Sequence_x, fill = CoclusteringRate)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_fill_steps(low = "white", high = "red", limit = c(0,1), space = "Lab", name = "Coclustering\nrate") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), legend.position = c(0.9, 0.9)) +
        ggplot2::coord_fixed() + ggplot2::xlab(NULL) +  ggplot2::ylab(NULL)
      filename <- paste(outputFolder, "/tilePlot_Cluster", clusIndex, ".", device, sep = "")
      do.call(ggplot2::ggsave, args = c(list(filename = filename, plot = plotObject, device = device), controlListGgsave))
    }
    tilePlotList <- lapply(largeClusters, FUN = produceOneClusTilePlot)
    returnObject <- tilePlotList
  } else {
    frameToPlot <- subset(frameToPlot, subset = CoclusteringRate > 1e-8)
    plotObject <-
      ggplot2::ggplot(data = frameToPlot, ggplot2::aes(x = Sequence_x, y = Sequence_y, fill = CoclusteringRate)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_steps(low = "white", high = "red", limit = c(0,1), space = "Lab", name = "Coclustering\nrate") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), legend.position = c(0.9, 0.9)) +
      ggplot2::coord_fixed() + ggplot2::xlab(NULL) +  ggplot2::ylab(NULL)
    freqTable <- table(frameToPlot$Sequence_x)
    singletons <- names(freqTable)[which(freqTable == 1)]
    frameToPlot <- subset(frameToPlot, subset = !(Sequence_x %in% singletons))
    seqsToKeep <- names(clusterIndices)[clusterIndices %in% largeClusters]
    plotObject <- plotObject + ggplot2::scale_x_discrete(limits = rev(seqsToKeep)) +  ggplot2::scale_y_discrete(limits = seqsToKeep) + ggplot2::theme(axis.text = ggplot2::element_text(colour = "blue", size = textScaleFactor))
    returnObject <- plotObject
    filename <- paste(outputFolder, "/clustersTilePlot.", device, sep = "")
    do.call(ggplot2::ggsave, args = list(filename = filename, plot = plotObject, device = device, controlListGgsave))
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

.convertToSiteNum <- function(numVec, totalNumberOfNucleo, removedSites) {
  sitesThatWereKept <- setdiff(1:totalNumberOfNucleo, unique(removedSites))
  sitesThatWereKept[numVec]
}
