#' ---
#' title:  plot_violine.R
#' author: Sebastian Mueller (sebm_at_posteo.de)
#' date:   2018-04-10
#' desc: Creating varios summary statistics 

# o   A delimited file containing information for each STAMP (before cut-off) on number of UMIs, number of Genes detected/captured and the number of NGS-reads
# o   A separate delimited file containing information for each STAMP (after cut-off) on number of UMIs, number of Genes detected/captured and the number of NGS-reads
# Example of the format could be as follows:
# STAMP id          Number of NGS-reads       Number of UMIs       Number of Genes Detected
# STAMP1           1000000                           50000                       6000
#' ---
### for debug
# If you wish to access the snakefile object first invoke snakemake and save the session automatically
# Since there are no debug flags to my knowledge, just uncomment the line below and run snakemake which
# creates an R object that can be loaded into a custom R session

setwd("~/workspace/code/dropSeqPipe/.test")
load("snakemake_create_summary_stats.rdata")
# load(file="R_image_create_summary_stats.rdata")

if (snakemake@config$DEBUG) {
 message("In debug mode: saving R objects to inspect later")
 save(snakemake, file="snakemake_create_summary_stats.rdata")
}

####/debug
library(dplyr) # Dataframe manipulation
library(Matrix) # Sparse matrices
library(stringr)
library(RColorBrewer)
library(devtools)
library(Seurat)
library(plotly)

# importing Seurat object

# creating environment so objects don't get overwritten upon loading
env_imported_r_objects <- new.env()
load(file = file.path(snakemake@input$R_objects), envir = env_imported_r_objects)
# attach
seuratobj <- env_imported_r_objects$seuratobj
meta.data <- seuratobj@meta.data

#median calculator
meta.data %>%
  group_by(orig.ident) %>%
  summarise(median_number_genes     = median(nGene),
            median_Counts_per_STAMP = median(nCounts),
            median_UMIs_per_STAMP   = median(nUMI),
            mean_UMI_per_Gene       = mean(umi.per.gene),
            mean_Ribo_frac          = mean(pct.Ribo),
            mean_Mito_frac          = mean(pct.mito),
            read_length             = median(read_length), # should be all the same anyway..
            expected_cells          = median(expected_cells), # should be all the same anyway..
            actual_number_barcodes  = n()) %>%
   write.csv(file.path(snakemake@output$stats)) #writes table for excel           n = n()) %>%
# highest, lowest count/UMI Stamp
# pre STAMP stats

# hist out goes into knee plots 
		# 'logs/{sample}_hist_out_cell.txt'
		# """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && BAMTagHistogram -m {params.memory}\
		# TAG=XC\
# https://hpc.nih.gov/apps/dropseq.html
# there is not hint in documentation on any filtering (only read quality), but the lowest count is >1 (stange!))

if (snakemake@config$DEBUG) {
 save.image(file="R_image_create_summary_stats.rdata")
}
