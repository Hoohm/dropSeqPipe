library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(viridis)

samples = snakemake@params$sample_names
batches = snakemake@params$batches
data = data.frame(matrix(nrow=length(samples), ncol=10))
colnames(data) = c('Sample','Batch','BC_dropped','UMI_dropped','BC_and_UMI_dropped','Trimmomatic_filtered','Unmapped','Uniquely_mapped','Multi_mapped','Total_reads')
data[,'Sample'] = samples
data[,'Batch'] = batches
for(i in 1:length(samples)){
  #Input files and variables
  STAR_output = read.table(snakemake@input$STAR_output[i], skip = 5, sep = '\t', fill = TRUE, stringsAsFactors = FALSE)
  sample = samples[i]
  #Read files
  BC_UMI_filtered_reads_left = read.table(snakemake@input$reads_left[i])
  trimmomatic_filtered_reads_left = read.table(snakemake@input$trimmomatic_filtered[i])/4
  BC_tagged_data = read.table(snakemake@input$BC_tagged[i], header = TRUE)
  
  data[i,'Trimmomatic_filtered'] = BC_UMI_filtered_reads_left - trimmomatic_filtered_reads_left
  
  #STAR output
  reads_in = as.numeric(STAR_output$V2[1])
  uniquely_mapped = as.numeric(STAR_output$V2[4])
  multi_mapped = as.numeric(STAR_output$V2[19])
  unmapped = reads_in - uniquely_mapped - multi_mapped
  
  data[i,'Total_reads'] = sum(BC_tagged_data$num_barcodes)
  ids = which(BC_tagged_data$num_failed_bases > snakemake@params$min_num_below_BC)
  num_reads_BC_tagged = sum(BC_tagged_data$num_barcodes[ids])
  UMI_tagged_data = read.table(snakemake@input$UMI_tagged[i], header = TRUE)
  ids = which(UMI_tagged_data$num_failed_bases > snakemake@params$min_num_below_UMI)
  num_reads_UMI_tagged = sum(UMI_tagged_data$num_barcodes[ids])

  intersection = num_reads_BC_tagged + num_reads_UMI_tagged + BC_UMI_filtered_reads_left - data[i,'Total_reads']
  total_tags_filtered = data[i,'Total_reads'] - BC_UMI_filtered_reads_left
  #Since we don't know how much of cell and UMI barcodes intersect in terms of tagging, we assume linear proportions and substract them.
  num_reads_BC_dropped = num_reads_BC_tagged - intersection
  num_reads_UMI_dropped = num_reads_UMI_tagged - intersection
  #num_reads_BC_dropped = round(num_reads_BC_tagged - (intersection/total_tags_filtered)*num_reads_BC_tagged)
  #num_reads_UMI_dropped = total_tags_filtered - num_reads_BC_dropped

  data[i,'BC_and_UMI_dropped'] = intersection
  data[i,'BC_dropped'] = num_reads_BC_dropped
  data[i,'UMI_dropped'] = num_reads_UMI_dropped
  data[i,'Uniquely_mapped'] = uniquely_mapped
  data[i,'Multi_mapped'] = multi_mapped
  data[i,'Unmapped'] = unmapped
}


data_long = melt(data, c('Sample', 'Batch'))
p1 = ggplot(subset(data_long, data_long$variable != 'Total_reads'), aes(x=Sample, y = value, fill = factor(variable)))
p1 = p1 + geom_histogram(stat = 'identity', binwidth = 1/length(samples))
p1 = p1 + theme(axis.text.x=element_text(angle = 90, hjust = 0))
p1 = p1 + labs(title=paste('Yield of all the reads for each category'), x='Samples', y='Number of reads')
p1 = p1 + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.position="none", plot.title = element_text(size = 20, face = "bold"))
p1 = p1 + facet_grid(~Batch, scales = "free")
p1 = p1 + scale_fill_viridis(discrete=TRUE, option='magma')
p1 = p1 + scale_y_continuous(labels = scales::scientific)

data_pct = data[,-c(1,2)]/data$Total_reads
data_pct = cbind(Sample=data[,'Sample'],Batch = data[,'Batch'],data_pct)

data_long = melt(data_pct, c('Sample', 'Batch'))

p2 = ggplot(subset(data_long, data_long$variable != 'Total_reads'), aes(x=Sample, y = value, fill = factor(variable, labels=c(
                                                                                      'Cell BC filtered',
                                                                                      'UMI BC filtered',
                                                                                      'Cell and UMI filtered',
                                                                                      'Trimmomatic filtered',
                                                                                      'Unmapped',
                                                                                      'Uniquely mapped',
                                                                                      'Multiply mapped')))) + labs(fill = "Filters")
p2 = p2 + geom_histogram(stat = 'identity', binwidth = 1/length(samples))
p2 = p2 + theme(axis.text.x=element_text(angle = 90, hjust = 0))
p2 = p2 + labs(x='Samples', y='Percentage of reads')
p2 = p2 + theme(legend.position="bottom")
p2 = p2 + facet_grid(~Batch, scales = "free")
p2 = p2  + scale_fill_viridis(discrete=TRUE, option='magma')
p2 = p2 + scale_y_continuous(labels = scales::percent)

gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)

pdf(file = snakemake@output$pdf, width = 16, height = 13)
grid::grid.newpage()
grid::grid.draw(rbind(gp1, gp2, size = "last"))
dev.off()