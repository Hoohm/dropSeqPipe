library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
samples = snakemake@params$sample_names
data = data.frame(matrix(nrow=length(samples), ncol=6))
colnames(data) = c('Sample','BC_dropped','UMI_dropped','BC_UMI_dropped','Not_dropped','Total_reads')
data[,'Sample'] = samples
#fastqc_summary = read.table(snakemake@input$fastqc_txt, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, sep='\t')
for(i in 1:length(samples)){
  sample = samples[i]
  not_dropped = read.table(snakemake@input$reads_left[i])
  data[i,'Not_dropped'] = not_dropped
  #data[i,'Total_reads'] = as.integer(unique(as.character(fastqc_summary[grep(sample, rownames(fastqc_summary)),'FastQC_total_sequences-1'])))
  BC_tagged_data = read.table(snakemake@input$BC_tagged[i], header = TRUE)
  data[i,'Total_reads'] = sum(BC_tagged_data$num_barcodes)
  #BC_tagged_data = BC_tagged_data[-1,] # delete 0 error line
  ids = which(snakemake@params$min_num_below_BC < BC_tagged_data$num_failed_bases)
  num_reads_BC_dropped = (data[i,'Total_reads'] - data[i,'Not_dropped']) - sum(BC_tagged_data$num_barcodes[ids])
  data[i,'BC_dropped'] = num_reads_BC_dropped
  UMI_tagged_data = read.table(snakemake@input$UMI_tagged[i], header = TRUE)
  #UMI_tagged_data = UMI_tagged_data[-1,]# delete 0 error line
  ids = which(snakemake@params$min_num_below_UMI < UMI_tagged_data$num_failed_bases)
  num_reads_UMI_dropped = (data[i,'Total_reads'] - data[i,'Not_dropped']) - sum(UMI_tagged_data$num_barcodes[ids])
  data[i,'UMI_dropped'] = num_reads_UMI_dropped
  data[i,'BC_UMI_dropped'] = (data[i,'Total_reads'] - data[i,'Not_dropped']) - sum(data[i,'UMI_dropped'] + data[i,'BC_dropped'])
}
data_long = melt(data, 'Sample')
# Keep the order of the barcodes using factor and levels.
p1 = ggplot(subset(data_long, data_long$variable != 'Total_reads'), aes(x=Sample, y = value, fill = variable))
p1 = p1 + geom_bar(stat = 'identity')
p1 = p1 + theme(axis.text.x=element_text(angle = 90, hjust = 0))
p1 = p1 + labs(title=paste('BC and UMI drop comparison.\nMin Cell BC quality =',snakemake@params$min_BC_quality,'& Max number bellow =',snakemake@params$min_num_below_BC,'\nMin UMI BC quality =',snakemake@params$min_UMI_quality,'& Max number bellow =',snakemake@params$min_num_below_UMI), x='Samples', y='Number of reads')
p1 = p1 + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p1 = p1 + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#10d09b"))
data$BC_dropped = data$BC_dropped/data$Total_reads
data$UMI_dropped = data$UMI_dropped/data$Total_reads
data$BC_UMI_dropped = data$BC_UMI_dropped/data$Total_reads
data_long_pct = melt(data,'Sample', measure.vars = c('BC_dropped', 'UMI_dropped', 'BC_UMI_dropped'))
p2 = ggplot(data_long_pct, aes(x=Sample, y = value, fill = variable))
p2 = p2 + geom_bar(stat = 'identity')
p2 = p2 + theme(axis.text.x=element_text(angle = 90, hjust = 0))
p2 = p2 + labs(x='Samples', y='Percentage of reads')
p2 = p2 + scale_y_continuous(labels = scales::percent)
p2 = p2 + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
# This allows to align the main plots so that we can relate both directly with the label from the bottom one.
gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)
pdf(file = snakemake@output[[1]], width = 8, height = 7)
grid::grid.newpage()
grid::grid.draw(rbind(gp1, gp2, size = "last"))
dev.off()