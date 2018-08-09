library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(viridis)

samples = snakemake@params$sample_names
batches = snakemake@params$batches
data = data.frame(matrix(nrow=length(samples), ncol=8))
colnames(data) = c('Sample','Batch','Cell_dropped','UMI_dropped','Cell_UMI_dropped','Trimmomatic_filtered','Not_dropped','Total_reads')
data[,'Sample'] = samples
data[,'Batch'] = batches

for(i in 1:length(samples)){
  sample = samples[i]
  #Read files
  Cell_UMI_filtered_reads_left = read.table(snakemake@input$reads_left[i])
  trimmomatic_filtered_reads_left = read.table(snakemake@input$trimmomatic_filtered[i])/4
  Cell_tagged_data = read.table(snakemake@input$Cell_tagged[i], header = TRUE)
  
  data[i,'Trimmomatic_filtered'] = Cell_UMI_filtered_reads_left - trimmomatic_filtered_reads_left
  data[i,'Total_reads'] = sum(Cell_tagged_data$num_barcodes)
  tagged_cell_ids = which(Cell_tagged_data$num_failed_bases > snakemake@params$min_num_below_Cell)
  num_reads_Cell_tagged = sum(Cell_tagged_data$num_barcodes[tagged_cell_ids])
  UMI_tagged_data = read.table(snakemake@input$UMI_tagged[i], header = TRUE)
  tagged_umi_ids = which(UMI_tagged_data$num_failed_bases > snakemake@params$min_num_below_UMI)
  num_reads_UMI_tagged = sum(UMI_tagged_data$num_barcodes[tagged_umi_ids])

  intersection = num_reads_Cell_tagged + num_reads_UMI_tagged + Cell_UMI_filtered_reads_left - data[i,'Total_reads']
  total_tags_filtered = data[i,'Total_reads'] - Cell_UMI_filtered_reads_left
  #Since we don't know how much of cell and UMI barcodes intersect in terms of tagging, we assume linear proportions and substract them.
  num_reads_Cell_dropped = num_reads_Cell_tagged - intersection
  num_reads_UMI_dropped = num_reads_UMI_tagged - intersection
  #num_reads_Cell_dropped = round(num_reads_Cell_tagged - (intersection/total_tags_filtered)*num_reads_Cell_tagged)
  #num_reads_UMI_dropped = total_tags_filtered - num_reads_Cell_dropped
  
  data[i,'Cell_UMI_dropped'] = intersection
  data[i,'Cell_dropped'] = num_reads_Cell_dropped
  data[i,'UMI_dropped'] = num_reads_UMI_dropped
  data[i,'Not_dropped'] = trimmomatic_filtered_reads_left
}

data_long = melt(data, c('Sample','Batch'))

# Keep the order of the barcodes using factor and levels.
p1 = ggplot(subset(data_long, data_long$variable != 'Total_reads'), aes(x=Sample, y = value, fill = variable))
p1 = p1 + geom_bar(stat = 'identity')
p1 = p1 + theme(axis.text.x=element_text(angle = 90, hjust = 0))
p1 = p1 + labs(title=paste('Cell and UMI drop comparison.\nMin Cell Cell quality =',snakemake@params$min_Cell_quality,'& Max number bellow =',snakemake@params$min_num_below_Cell,'\nMin UMI Cell quality =',snakemake@params$min_UMI_quality,'& Max number bellow =',snakemake@params$min_num_below_UMI), x='Samples', y='Number of reads')
p1 = p1 + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p1 = p1 + facet_grid(~Batch, scales = "free")
p1 = p1 + scale_fill_viridis(discrete=TRUE, option='magma')
p1 = p1 + scale_y_continuous(labels = scales::scientific)
#p1 = p1 + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#10d09b"))


data[,-c(1,2)] = data[,-c(1,2)]/data$Total_reads

data_long_pct = melt(data,c('Sample', 'Batch'), measure.vars = c('Cell_dropped', 'UMI_dropped','Cell_UMI_dropped', 'Trimmomatic_filtered', 'Not_dropped'))
p2 = ggplot(data_long_pct, aes(x=Sample, y = value, fill = variable))
p2 = p2 + geom_bar(stat = 'identity')
p2 = p2 + theme(axis.text.x=element_text(angle = 90, hjust = 0))
p2 = p2 + labs(x='Samples', y='Percentage of reads')
p2 = p2 + scale_y_continuous(labels = scales::percent)
p2 = p2 + facet_grid(~Batch, scales = "free")
p2 = p2 + scale_fill_viridis(discrete=TRUE, option='magma')
# This allows to align the main plots so that we can relate both directly with the label from the bottom one.
gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)
pdf(file = snakemake@output$pdf, width = 8, height = 7)
grid::grid.newpage()
grid::grid.draw(rbind(gp1, gp2, size = "last"))
dev.off()
