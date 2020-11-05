sink(snakemake@log[["stdout"]])

debug_flag <- FALSE
if (snakemake@config$DEBUG) {
  debug_flag <- TRUE
  message("In debug mode: saving R objects to inspect later")
  path_debug <- file.path(snakemake@config$LOCAL$results, "debug")
  dir.create(path_debug, showWarnings = FALSE)
  save(snakemake, file = file.path(path_debug, "create_seurat_object_snakemake.rdata"))
}


options(warn = -1)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE) # Dataframe manipulation
library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(Seurat, quietly = TRUE, warn.conflicts = FALSE)
library(DropletUtils, quietly = TRUE, warn.conflicts = FALSE)

design <- read.csv(snakemake@params$design)

sce_list <- list()
sce_list_filtered <- list()
n=1
for (sample_path in snakemake@input$umi_mtx) {
  # sample_path = snakemake@input$umi_mtx[n]
  sample_name <- snakemake@params$samples[n]
  print(sample_name)
  print(sample_path)
  if (!grepl(pattern = sample_name, x = sample_path)) {
    stop("sample order is not proper. Exiting")
  }
  # temp_data = Read10X(dirname(sample_path))
  # using DropletUtils to read in raw data to use emptyDrops for filtering
  sce <- DropletUtils::read10xCounts(dirname(sample_path), version = "3",
                                     type = "sparse", col.names = TRUE)
  sce_list[sample_name] <- sce
  barcoderanks <- barcodeRanks(counts(sce))

  png(paste0(sample_name, ".barcoderank.png"))
  plot(barcoderanks$rank, barcoderanks$total, log = "xy",
       xlab = "Rank", ylab = "Total")
  o <- order(barcoderanks$rank)
  lines(barcoderanks$rank[o], barcoderanks$fitted[o], col = "red")
  abline(h = metadata(barcoderanks)$knee, col = "dodgerblue", lty = 2)
  abline(h = metadata(barcoderanks)$inflection, col = "forestgreen", lty = 2)
  legend("bottomleft",
    lty = 2, col = c("dodgerblue", "forestgreen"),
    legend = c("knee", "inflection")
  )
  dev.off()

  # https://bioconductor.org/packages/devel/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html#detecting-empty-droplets
  set.seed(42)
  # e.out <- emptyDrops(counts(sce), lower = 10, retain = 1000)
  # TODO: offering option to use knee-plot threshold or expected-cells instead of emptyDrops
  e.out <- emptyDrops(counts(sce))
  print(e.out)

  is.cell <- e.out$FDR <= 0.01
  sum(is.cell, na.rm = TRUE)

  table(Limited = e.out$Limited, Significant = is.cell)

  is.cell2 <- is.cell
  is.cell2[is.na(is.cell2)] <- FALSE
  sce_filtered <- sce[, is.cell2]

  png(paste0(sample_name, ".logProb_vs_UMIrank.png"))
  plot(e.out$Total, -e.out$LogProb,
    col = ifelse(is.cell, "red", "black"),
    xlab = "Total UMI count", ylab = "-Log Probability"
  )
  dev.off()

  if (ncol(sce_filtered) > 0) {
    # create SingleCellExperiment object
    sce_list_filtered[sample_name] <- sce_filtered
    # create Seurat object
    seuratobj <- CreateSeuratObject(counts = assay(sce_filtered, "counts"), project = sample_name)
  } else {
    message("no cells left after filtering step for ", sample_name)
  }

  n <- n + 1

  expdat <- assay(sce_filtered, "counts")
}

# Merge all if more than 1 sample
if (length(sce_list_filtered) > 1) {
  sce_all_filtered <- Reduce(SingleCellExperiment::cbind, sce_list_filtered[-1], sce_list_filtered[[1]])
  sce_all <- Reduce(SingleCellExperiment::cbind, sce_list[-1], sce_list[[1]])
} else {
  sce_all_filtered <- sce_list_filtered[[1]]
  sce_all <- sce_list[[1]]
}

# TODO creating Seurat object and adding rna_metrics infor
    # this still needs tweaking since barcodes are non-unique and orig.ident is wrong
    # seurat_obj = merge(seurat_obj, temp_object)


  # rna_metrics_data <- read.csv(file = snakemake@input$rna_metrics[n], header = T,
  #                    stringsAsFactors = F, skip = 6, sep = "\t")
  # meta_data = rna_metrics_data %>% filter(SAMPLE %in% colnames(temp_data))

  # sub_design = design %>% filter(samples==sample_name)
  # for (column in colnames(sub_design)){
  #   if (column != 'samples'){
  #     meta_data[,column] = sub_design[,column]
  #   }
  # }
  # rownames(meta_data) = rna_metrics_data$SAMPLE
  # colnames(temp_data) = paste0(colnames(temp_data), '-', n)
  # rownames(meta_data) = paste0(rownames(meta_data), '-', n)

  # temp_object = CreateSeuratObject(counts = temp_data, project = "sample_name")
  # if (n==1){
  #   seurat_obj = temp_object
  # } else{
  #   seurat_obj = merge(seurat_obj, temp_object)
  # }
  # n=n+1

# seurat_obj@misc[['config']] = snakemake@config
seurat_obj <- "placeholder"
# create dir if it doesn't exist
dir.create(dirname(snakemake@output$seurat_object))
saveRDS(seurat_obj, file = snakemake@output$seurat_object)
# saving sce object instead of Seurat for the time being
saveRDS(sce_all_filtered, file = snakemake@output$SCE_object)

sessionInfo()
sink()
