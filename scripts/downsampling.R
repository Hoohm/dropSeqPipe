#------------------------------------ for debugging:
# For debugging add the following line in config.yaml (without the #)
# DEBUG: True
# This will create R objects in the debug directory containing the snakemake
# object R object that can be loaded into a custom R session as below:
debug_flag <- FALSE
if (snakemake@config$DEBUG) {
  debug_flag <- TRUE
  message("In debug mode: saving R objects to inspect later")
  path_debug <- file.path(snakemake@config$LOCAL$results, "debug")
  dir.create(path_debug)
  save(snakemake, file = file.path(path_debug, "subsampling_snakemake.rdata"))
}
#------------------------------------ debugging


library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)


results = NULL
filtered_data = 'data/filtered/run13_new_filtered.Rds'
filtered = readRDS(file = filtered_data)



bc_umi = read.table(snakemake@input$bc_umi), sep = '\t', header = TRUE)



n_cells = 
k_reads_per_cell = c(seq(10000,90000, 10000), seq(100000, 1000000, 100000))

for (sample in samples){
for (i in k_reads_per_cell){
  #for (j in 1:10){
  #sample_data = subset[,ids]
  k_reads = i*n_cells
  total_reads = sum(subset$Num_Obs)
  if(k_reads < total_reads){
    print(i)
    #sample_data = as.data.frame(sample_data)
  #sample_data$genes = rownames(sample_data)
  
  #long = sample_data %>% gather(key=cell, value=expr, -genes)
  
  downsampled = subset %>% sample_n(size = k_reads, replace = TRUE, weight = Num_Obs) %>% group_by(Cell.Barcode, Gene) %>% distinct()  %>% summarise(umi=n())%>% spread(Cell.Barcode, umi, fill = 0) %>% ungroup()  %>% select(-Gene)
  counts_per_genes = rowSums(downsampled)
  n_genes = sum(counts_per_genes!=0)
  n_umis = sum(counts_per_genes)
  umi_per_cell = colSums(downsampled)
  per_cell_genes = apply(downsampled, 2, function(x) {sum(x!=0)})
  mean_umi_per_cell = mean(umi_per_cell[which(umi_per_cell!=0)])#apply(downsampled, 2, function(x) {mean(x[which(x!=0)])})
  median_umi_per_cell = median(umi_per_cell[which(umi_per_cell!=0)])#apply(downsampled, 2, function(x) {median(x[which(x!=0)])})
  
  mean_genes = mean(per_cell_genes, na.rm=TRUE)
  median_genes = median(per_cell_genes, na.rm=TRUE)
  
  
  
  
  results = rbind(results, data.frame(i, n_genes, n_umis, mean_genes, median_genes, mean_umi_per_cell, median_umi_per_cell, sample))
  #}
  } 
}
  
  
  
}
print('done')
}

colnames(results) = c('k_filtered_reads_per_cell','total_genes','total_umis', 'mean_genes_per_cell', 'median_genes_per_cell', 'mean_umis_per_cell','median_umis_per_cell','sample')

p1 = ggplot(results, aes(x=k_filtered_reads_per_cell, y=median_genes_per_cell, color=sample)) +geom_line() + theme(legend.position = 'none') +xlim(10000,90000)
p2 = ggplot(results, aes(x=k_filtered_reads_per_cell, y=mean_genes_per_cell, color=sample)) +geom_line() + theme(legend.position = 'none') +xlim(10000,90000)

p3 = ggplot(results, aes(x=k_filtered_reads_per_cell, y=median_umis_per_cell, color=sample)) +geom_line() + theme(legend.position = 'none') +xlim(10000,90000)
p4 = ggplot(results, aes(x=k_filtered_reads_per_cell, y=mean_umis_per_cell, color=sample)) +geom_line() + theme(legend.position = 'none') +xlim(10000,90000)
p5 = ggplot(results, aes(x=k_filtered_reads_per_cell, y=total_genes, color=sample)) +geom_line() + theme(legend.position = 'none') +xlim(10000,90000)
p6 = ggplot(results, aes(x=k_filtered_reads_per_cell, y=total_umis, color=sample)) +geom_line() +xlim(10000,90000)
p1+p2+p3+p4+p5+p6 +plot_layout(ncol=2)

if (debug_flag) {
  save.image(file = file.path(path_debug, "plot_yield_workspace.rdata"))
}
