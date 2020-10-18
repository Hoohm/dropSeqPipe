
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

design = read.csv(snakemake@params$design)

n=1
for (sample_path in snakemake@input$umi_mtx){
  sample_name = snakemake@params$samples[n]
  print(sample_name)
  print(sample_path)
  if (!grepl(pattern = sample_name, x = sample_path)){
    quit('sample order is not proper. Exiting')
  }
  temp_data = Read10X(dirname(sample_path))
  
  rna_metrics_data <- read.csv(file = snakemake@input$rna_metrics[n], header = T,
                     stringsAsFactors = F, skip = 6, sep = "\t")
  meta_data = rna_metrics_data %>% filter(SAMPLE %in% colnames(temp_data))
  
  sub_design = design %>% filter(samples==sample_name)
  for (column in colnames(sub_design)){
    if (column != 'samples'){
      meta_data[,column] = sub_design[,column]
    }
  }
  rownames(meta_data) = rna_metrics_data$SAMPLE
  colnames(temp_data) = paste0(colnames(temp_data), '-', n)
  rownames(meta_data) = paste0(rownames(meta_data), '-', n)
  
  temp_object = CreateSeuratObject(counts = temp_data, project = sample_name, meta.data = meta_data)
  if (n==1){
    seurat_obj = temp_object  
  } else{
    seurat_obj = merge(seurat_obj, temp_object)  
  }
  n=n+1
  
}


seurat_obj@misc[['config']] = snakemake@config

save(seurat_obj, file = snakemake@output$seurat_object)
