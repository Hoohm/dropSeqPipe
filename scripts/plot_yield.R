library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(viridis)
library(stringr)

samples = snakemake@params$sample_names
batches = snakemake@params$batches
data = data.frame(matrix(nrow=length(samples), ncol=7))
colnames(data) = c('Sample','Batch','Cutadapt filtered','Unmapped','Multi mapped','Uniquely mapped','Total reads')
data[,'Sample'] = samples
data[,'Batch'] = batches
for(i in 1:length(samples)){
  #Input files and variables
  STAR_output = read.table(snakemake@input$STAR_output[i], skip = 5, sep = '\t', fill = TRUE, stringsAsFactors = FALSE)
  sample = samples[i]
  #Read files
  bbmap_log = read.table(snakemake@input$repaired[i], sep=':', header=FALSE, skip=6, row.names=1, nrows=4)
  reads_after_filtering = as.numeric(str_match(bbmap_log['Pairs',], pattern = "\t([0-9]{1,20}) reads.*")[,2])/2
  R1_filtered = read.table(snakemake@input$R1_filtered[i], header = FALSE, skip=8, sep=':', nrows=7, row.names=1)
  total_reads = as.numeric(str_replace_all(R1_filtered['Total reads processed',], pattern = (' |,'), ""))
  
  # R2_filtered = read.table(snakemake@input$R2_filtered[i], header = FALSE, skip=8, sep=':', nrows=7, row.names=1)
  
  # R1_adapters = as.numeric(str_remove_all(str_match(R1_filtered['Reads with adapters',], pattern = "(.*) \\(")[,2], pattern = (' |,')))
  # R1_too_short = as.numeric(str_remove_all(str_match(R1_filtered['Reads that were too short',], pattern = "(.*) \\(")[,2], pattern = (' |,')))
  # R1_passed = as.numeric(str_remove_all(str_match(R1_filtered['Reads written (passing filters)',], pattern = "(.*) \\(")[,2], pattern = (' |,')))
  # R1_filtered = total_reads - R1_passed
  
  # R2_adapters = as.numeric(str_remove_all(str_match(R2_filtered['Reads with adapters',], pattern = "(.*) \\(")[,2], pattern = (' |,')))
  # R2_too_short = as.numeric(str_remove_all(str_match(R2_filtered['Reads that were too short',], pattern = "(.*) \\(")[,2], pattern = (' |,')))
  # R2_passed = as.numeric(str_remove_all(str_match(R2_filtered['Reads written (passing filters)',], pattern = "(.*) \\(")[,2], pattern = (' |,')))
  # R2_filtered = total_reads - R2_passed


  data[which(data$Sample == sample),'Cutadapt filtered'] = total_reads - reads_after_filtering
  data[which(data$Sample == sample),'Total reads'] = total_reads

  
  

  
  #STAR output
  reads_in = as.numeric(STAR_output$V2[1])
  uniquely_mapped = as.numeric(STAR_output$V2[4])
  multi_mapped = as.numeric(STAR_output$V2[19])
  unmapped = reads_in - uniquely_mapped - multi_mapped
  
  data[which(data$Sample == sample),'Uniquely mapped'] = uniquely_mapped
  data[which(data$Sample == sample),'Multi mapped'] = multi_mapped
  data[which(data$Sample == sample),'Unmapped'] = unmapped
  
}

data_long = melt(data, c('Sample', 'Batch'))

p1 = ggplot(subset(data_long, data_long$variable != 'Total reads'), aes(x=Sample, y = value, fill = factor(variable)))
p1 = p1 + geom_histogram(stat = 'identity', binwidth = 1/length(samples))
p1 = p1 + theme(axis.text.x=element_text(angle = 90, hjust = 0))
p1 = p1 + labs(title=paste('Yield of all the reads for each category'), x='Samples', y='Number of reads')
p1 = p1 + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.position="none", plot.title = element_text(size = 20, face = "bold"))
p1 = p1 + facet_grid(~Batch, scales = "free")
p1 = p1 + scale_fill_viridis(discrete=TRUE, option='magma')
p1 = p1 + scale_y_continuous(labels = scales::scientific)

data_pct = data[,-c(1,2)]/data[,'Total reads']
data_pct = cbind(Sample=data[,'Sample'],Batch = data[,'Batch'],data_pct)

data_long = melt(data_pct, c('Sample', 'Batch'))

p2 = ggplot(subset(data_long, data_long$variable != 'Total reads'), aes(x=Sample, y = value, fill = variable)) + labs(fill = "Filters")
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
