#' ---
#' title:  plot_violine.R
#' author: Sebastian Mueller (sebm_at_posteo.de)
#' date:   2018-04-10
#' ---
### for debug
# If you wish to access the snakefile object first invoke snakemake and save the session automatically
# Since there are no debug flags to my knowledge, just uncomment the line below and run snakemake which
# creates an R object that can be loaded into a custom R session
# save.image(file="R_workspace_debug.rdata")
# load("R_workspace_debug.rdata")
####/debug
library(plyr)
library(dplyr) # Dataframe manipulation
library(Matrix) # Sparse matrices
library(stringr)
library(RColorBrewer)
#library(matrixStats)
# library(Hmisc) # for cut2 function
library(devtools)
library(Seurat)
library(plotly)

# rule map in Snakefile
# rule map:
#     input:
#         'plots/violinplots_comparison_UMI.pdf',
#         ...

# importing UMI
# umi_matrix             <- read.csv(file.path(path,'summary/umi_expression_matrix.tsv'), sep='\t', header = TRUE, row.names = 1, check.names = FALSE) %>%
# importing counts ( summary/counts_expression_matrix.tsv )
count_matrix <- read.csv(snakemake@input$counts, sep = "\t",
                         header = TRUE, row.names = 1,
                         check.names = FALSE)
# importing UMIs ( summary/umi_expression_matrix.tsv )
umi_matrix   <- read.csv(snakemake@input$UMIs,
                         sep = "\t",
                         header = TRUE,
                         row.names = 1,
                         check.names = FALSE)
design       <- read.csv(snakemake@input$design, stringsAsFactors = TRUE,
                         header = TRUE,
                         row.names = NULL)

metaData <- data.frame(cellNames = colnames(umi_matrix)) %>%
  mutate(samples = factor(str_replace(cellNames,"_[^_]*$",""))) %>%
  mutate(barcode = factor(str_replace(cellNames,".+_",""))) %>%
  left_join(design, by = "samples")
rownames(metaData) <- metaData$cellNames

# possible to set is.expr = -1 to avoid filtering whilst creating
# myumi <- CreateSeuratObject(raw.data = umi_matrix, meta.data = metaData, is.expr = -1)
myumi <- CreateSeuratObject(raw.data = umi_matrix, meta.data = metaData)
myumi <- SetAllIdent(object = myumi, id = "samples")
# relabel cell idenity (https://github.com/satijalab/seurat/issues/380)
myumi@meta.data$orig.ident <- myumi@meta.data$samples

mycount <- CreateSeuratObject(raw.data = count_matrix, meta.data = metaData)
mycount <- SetAllIdent(object = mycount, id = "samples")
mycount@meta.data$orig.ident <- mycount@meta.data$samples
# turn off filtering

# note, the @meta.data slot contains usefull summary stuff
# head(mycount@meta.data,2)
#                              nGene nUMI expected_cells read_length      barcode
# dropseqLib1_ACTAACATTATT    15   33            400         100 ACTAACATTATT
# dropseqLib1_GAGTCTGAGGCG     5    9            400         100 GAGTCTGAGGCG
#                                       origin      origin
# dropseqLib1_ACTAACATTATT dropseqLib1 dropseqLib1
# dropseqLib1_GAGTCTGAGGCG dropseqLib1 dropseqLib1
meta.data         <- myumi@meta.data
meta.data$nCounts <- mycount@meta.data$nUMI


# mytheme <- theme_bw(base_size = 9) +
mytheme <- theme_bw() +
  theme(legend.position = "right",
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 300, hjust = 0))
theme_set(mytheme)

# predefined ggplot layers for subsequent plots
gglayers <- list(
  geom_smooth(),
  geom_point(size = .5),
  scale_y_continuous(labels = scales::unit_format(unit = "", scale = 1e-3, digits = 2),
                     breaks = scales::pretty_breaks(n = 8)),
  scale_x_continuous(labels = scales::unit_format(unit = "", scale = 1e-3, digits = 2),
                     breaks = scales::pretty_breaks(n = 8))
)

gg <- ggplot(meta.data, aes(x = nUMI, y=nCounts, color=orig.ident)) +
  #   coord_trans(y="log10",x = "log10") +
  gglayers +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "UMI counts vs raw Counts",
     subtitle = "Number of UMIs and raw Counts for each Bead",
     x = "Number of UMIs per Bead [k]",
     y = "Number of Counts per Bead [k]")

htmlwidgets::saveWidget(ggplotly(gg), file.path(getwd(),snakemake@output$html_umivscounts))
ggsave(gg, file = file.path(getwd(), snakemake@output$pdf_umivscounts), width=12,height=7)

# how about unaligned reads/UMI?
# Note(Seb): raw.data is actually filtered data i.e. nr of genes likely to be smaller than input data!
mito.gene.names  <- grep("^mt-", rownames(myumi@raw.data), value=TRUE)
# mito.gene.names2  <- subset(mito.gene.names, mito.gene.names %in% rownames(myumi@raw.data))
sribo.gene.names <- grep("^Rps", rownames(myumi@raw.data), value=TRUE)
lribo.gene.names <- grep("^Rpl", rownames(myumi@raw.data), value=TRUE)

col.total.umi            <- Matrix::colSums(myumi@raw.data)
col.total.count          <- Matrix::colSums(mycount@raw.data)
meta.data$col.total.umi   <- col.total.umi
meta.data$col.total.count <- col.total.count
# mito.percent.counts <- mito.percent.counts
# sribo.pct <- sribo.pct
# lribo.pct <- lribo.pct
# ribo_tot <- ribo_tot

myumi.top_50   <- apply(myumi@raw.data, 2, function(x) sum(x[order(x, decreasing = TRUE)][1:50])/sum(x))
mycount.top_50 <- apply(mycount@raw.data, 2, function(x) sum(x[order(x, decreasing = TRUE)][1:50])/sum(x))

myumi <- AddMetaData(myumi, Matrix::colSums(myumi@raw.data[sribo.gene.names, ])/col.total.umi, "pct.sribo")
myumi <- AddMetaData(myumi, Matrix::colSums(myumi@raw.data[lribo.gene.names, ])/col.total.umi, "pct.lribo")
myumi <- AddMetaData(myumi, Matrix::colSums(myumi@raw.data[unique(c(sribo.gene.names, lribo.gene.names)), ])/col.total.umi, "pct.Ribo")
myumi <- AddMetaData(myumi, Matrix::colSums(myumi@raw.data[mito.gene.names, ])/col.total.umi, "pct.mito")
myumi <- AddMetaData(myumi, myumi.top_50, "top50")
tmp <- myumi@meta.data$nUMI/myumi@meta.data$nGene
names(tmp) <- rownames(myumi@meta.data)
myumi <- AddMetaData(myumi, tmp, "umi.per.gene")

# tmp=(myumi@meta.data[,"nUMI",drop=F]/myumi@meta.data$nGene)
mycount <- AddMetaData(mycount, Matrix::colSums(mycount@raw.data[sribo.gene.names, ])/col.total.count, "pct.sribo")
mycount <- AddMetaData(mycount, Matrix::colSums(mycount@raw.data[lribo.gene.names, ])/col.total.count, "pct.lribo")
mycount <- AddMetaData(mycount, Matrix::colSums(mycount@raw.data[unique(c(sribo.gene.names, lribo.gene.names)), ])/col.total.count, "pct.Ribo")
mycount <- AddMetaData(mycount, Matrix::colSums(mycount@raw.data[mito.gene.names, ])/col.total.count, "pct.mito")
mycount <- AddMetaData(mycount, mycount.top_50, "top50")
tmp <- mycount@meta.data$nUMI/mycount@meta.data$nGene
names(tmp) <- rownames(mycount@meta.data)
mycount <- AddMetaData(mycount, tmp, "umi.per.gene")
# mycount@meta.data$count.per.gene <- mycount@meta.data$nUMI/mycount@meta.data$nGene


gg <- VlnPlot(myumi,
              c("nUMI", "nGene", "top50", "umi.per.gene", "pct.Ribo", "pct.mito"),
              x.lab.rot = TRUE, do.return = TRUE)
# ggsave(gg,file=file.path("violinplots_comparison_UMI.pdf"),width=18,height=18)
ggsave(gg, file  = snakemake@output$pdf_violine, width = 18, height = 18)
# gg <- VlnPlot(mycount,c("nUMI", "nGene", "top50", "count.per.gene","pct.Ribo", "pct.mito"), x.lab.rot = TRUE, do.return = TRUE)
# ggsave(gg,file=file.path("violinplots_comparison_count.pdf"),width=18,height=18)

# gg <- GenePlot(object = myumi, gene1 = "nUMI", gene2 = "nGene")
# ggsave(gg,file=file.path("violinplots_comparison.pdf"),width=18,height=18)


gg <- ggplot(myumi@meta.data, aes(x = nUMI, y = nGene, color=orig.ident)) +
  gglayers +
  labs(title = "Genes (pooled mouse and human set) vs UMIs for each bead",
       x = "Number of UMIs per Bead [k]",
       y = "Number of Genes per Bead [k]")

htmlwidgets::saveWidget(ggplotly(gg),
                        file.path(getwd(), snakemake@output$html_umi_vs_gene))
# ggsave(gg, file = file.path(getwd(), snakemake@output$pdf_umi_vs_gene),
#        width = 12, height = 7)



################################################################################
## same for Counts instead UMIs (using mycount object)
gg <- ggplot(mycount@meta.data, aes(x = nUMI, y = nGene, color=orig.ident)) +
  gglayers +
  labs(title = "Genes (pooled mouse and human set) vs Counts for each bead",
       x = "Number of Counts per Bead [k]",
       y = "Number of Genes per Bead [k]")

htmlwidgets::saveWidget(ggplotly(gg),
                        file.path(getwd(), snakemake@output$html_count_vs_gene))

# ggsave(gg, file = file.path(getwd(), snakemake@output$pdf_count_vs_gene),
#        width = 12, height = 7)


save(snakemake, myumi, mycount,
     file=file.path(getwd(), snakemake@output$R_objects))
