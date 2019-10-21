# Functions used to plot the species plot for drop-seq mixed protocol
# Authors: James Nemesh, Roelli Patrick, Sebastian Y Mueller

debug_flag <- FALSE
if (snakemake@config$DEBUG) {
  debug_flag <- TRUE
  message("In debug mode: saving R objects to inspect later")
  path_debug <- file.path(snakemake@config$LOCAL$results, "debug")
  dir.create(path_debug, showWarnings = FALSE)
  save(snakemake, file = file.path(
    path_debug,
    paste0("plot_species_plot_snakemake_", attr(snakemake, "wildcard")$sample, ".rdata")
  ))
}

#### /debug

categorizeCellsUsingKneeKnownNumCellsPaper <- function(
                                                   digitalExpressionFileO1,
                                                   digitalExpressionFileO2,
                                                   organismOne,
                                                   organismTwo,
                                                   pureRatio = 0.2,
                                                   numCells,
                                                   numBeads,
                                                   point.cex = 1.5,
                                                   xlim_range = NULL,
                                                   ylim_range = NULL,
                                                   category = "transcripts") {
  dfFull <- getNumTranscriptsPerCellBarcodeByOrganismPair(
                                                    digitalExpressionFileO1,
                                                    digitalExpressionFileO2,
                                                    organismOne,
                                                    organismTwo,
                                                    category)
  dfFull <- dfFull[order(dfFull$total, decreasing = T), ]
  dfFull$ratio_one <- dfFull[, 2] / dfFull[, 4]
  dfFull <- head(dfFull, n = numBeads)
  df <- head(dfFull, n = numCells)

  dfNoCall <- dfFull[-1:-numCells, ]
  if (dim(dfNoCall)[1] > 0) {
    dfNoCall$organism <- "No Call"
  }

  df$organism <- "Mixed"

  idx <- which(df$ratio_one >= (1 - pureRatio))
  # checks if the species is actually assigned at all
  if (length(idx) > 0) {
    df[idx, ]$organism <- organismOne
  }
  idx <- which(df$ratio_one <= (pureRatio))
  if (length(idx) > 0) {
    df[idx, ]$organism <- organismTwo
  }

  result <- rbind(df, dfNoCall)

  maxRange <- max(result[, 2], result[, 3])

  dforganismOne <- result[result$organism == organismOne, ]
  dforganismTwo <- result[result$organism == organismTwo, ]
  dfMixed <- result[result$organism == "Mixed", ]
  dfNoCall <- result[result$organism == "No Call", ]

  if (is.null(xlim_range)) {
    xlim_range <- c(0, maxRange)
  }

  if (is.null(ylim_range)) {
    ylim_range <- c(0, maxRange)
  }
  colors <- c("blue", "red", "purple", "grey")
  plot(dforganismOne[, 2], dforganismOne[, 3], col = colors[1],
       pch = 16, xlim = xlim_range, ylim = ylim_range,
       xlab = paste(organismOne, category),
       ylab = paste(organismTwo, category),
       cex = point.cex)
  points(dforganismTwo[, 2], dforganismTwo[, 3],
         col = colors[2], pch = 16, cex = point.cex)
  points(dfMixed[, 2], dfMixed[, 3],
         col = colors[3], pch = 16, cex = point.cex)
  points(dfNoCall[, 2], dfNoCall[, 3],
         col = colors[4], pch = 16, cex = point.cex)
  l <- c(paste(organismOne, dim(dforganismOne)[1]),
         paste(organismTwo, dim(dforganismTwo)[1]),
         paste("Mixed", dim(dfMixed)[1]),
         paste("No Call", dim(dfNoCall)[1]))
  legend("topright", legend = l, fill = colors)
  title(paste("Species plot based on", category))
  return(df)
}

getNumTranscriptsPerCellBarcodeByOrganismPair <- function(
                                                      digitalExpressionFileO1,
                                                      digitalExpressionFileO2,
                                                      organismOne,
                                                      organismTwo,
                                                      category) {
  if (is.null(organismOne) || is.null(organismTwo)) {
    return(NULL)
  }

  o1 <- getGenesAndTranscriptsPerCellBarcode(digitalExpressionFileO1)
  o2 <- getGenesAndTranscriptsPerCellBarcode(digitalExpressionFileO2)

  commonBC <- union(o1$cellBC, o2$cellBC)
  o1p <- o1[match(commonBC, o1$cellBC), ]
  o2p <- o2[match(commonBC, o2$cellBC), ]
  if (category == "genes") {
    df <- data.frame(tag = commonBC, o1Count = o1p$numGenes,
                     o2Count = o2p$numGenes, stringsAsFactors = F)
  }
  else {
    df <- data.frame(tag = commonBC, o1Count = o1p$numTranscripts,
                     o2Count = o2p$numTranscripts, stringsAsFactors = F)
  }

  idx1 <- which(is.na(df$o1Count))
  idx2 <- which(is.na(df$o2Count))
  if (length(idx1) > 0) df[idx1, ]$o1Count <- 0
  if (length(idx2) > 0) df[idx2, ]$o2Count <- 0

  df$total <- apply(df[, 2:3], 1, sum, na.rm = T)
  df <- df[order(df$total, decreasing = T), ]
  colnames(df)[2] <- organismOne
  colnames(df)[3] <- organismTwo
  return(df)
}


getGenesAndTranscriptsPerCellBarcode <- function(digitalExpressionFile) {
  a <- read.table(digitalExpressionFile, header = T, stringsAsFactors = F)
  colnames(a) <- c("cellBC", "numGenicReads", "numTranscripts", "numGenes")
  return(a)
}

digitalExpressionFileO1 <- snakemake@input[[1]][1]
digitalExpressionFileO2 <- snakemake@input[[2]][1]

num_cells <- snakemake@params$expected_cells

organismOne <- names(snakemake@config$META$species)[1]
organismTwo <- names(snakemake@config$META$species)[2]

par(mar = c(5, 4, 4, 2) + 0.5)

pdf(snakemake@output$genes_pdf, height = 8, width = 8)
df_temp <- categorizeCellsUsingKneeKnownNumCellsPaper(
  digitalExpressionFileO1,
  digitalExpressionFileO2,
  organismOne = organismOne,
  organismTwo = organismTwo,
  pureRatio = snakemake@config$META$ratio,
  numCells = num_cells,
  numBeads = num_cells * 2,
  point.cex = 1,
  category = "genes"
)
dev.off()



pdf(snakemake@output$transcripts_pdf, height = 8, width = 8)
df <- categorizeCellsUsingKneeKnownNumCellsPaper(
  digitalExpressionFileO1,
  digitalExpressionFileO2,
  organismOne = organismOne,
  organismTwo = organismTwo,
  pureRatio = snakemake@config$META$ratio,
  numCells = num_cells,
  numBeads = num_cells * 2,
  point.cex = 1,
  category = "transcripts"
)
dev.off()
organism1 <- subset(df, df$organism == organismOne)
organism2 <- subset(df, df$organism == organismTwo)

write.table(organism1$tag, snakemake@output$barcodes_species[1],
            row.names = F, col.names = F, quote = F)
write.table(organism2$tag, snakemake@output$barcodes_species[2],
            row.names = F, col.names = F, quote = F)

# save.image(paste0(snakemake@output$genes_pdf,".rdata"))
if (debug_flag) {
  save.image(file = file.path(
    path_debug,
    paste0(
      "plot_species_plot_workspace_",
      attr(snakemake, "wildcard")$sample, ".rdata"
    )
  ))
}
