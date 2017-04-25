library(ggplot2)
library(reshape2)
library(ggseqlogo)
library(yaml)
library(gridExtra)
library(grid)
library(ggpubr)
args = commandArgs(TRUE)
path = args[1]
config_file_data = yaml.load_file(paste0(path,'config.yaml'))
samples = names(config_file_data$Samples)


# Create the RNAmetrics plot
plotRNAMetrics = function(file_path, sample_name, top_barcodes, path){
  data = read.table(file = file_path, header=T, stringsAsFactors=F)
  data = data[order(data$PF_ALIGNED_BASES, decreasing = T),]
  # TODO: change to named columns
  data_pct = data[,c(25,11,12,13,14,15)]
  data = data[,c(25,3,4,5,6,7)]
  
  data_long = melt(data, id.var = "READ_GROUP")
  # Keep the order of the barcodes using factor and levels.
  my_sequences = levels(factor(data_long$READ_GROUP))
  data_long$READ_GROUP <- factor(data_long$READ_GROUP, levels = data_long$READ_GROUP)
  p1 = ggplot(data_long, aes(x=READ_GROUP, y = value, fill = variable)) + geom_bar(stat = 'identity') + theme(axis.text.x=element_text(angle = 90, hjust = 0)) + labs(title=paste('Top',top_barcodes, 'barcodes for', sample_name), x='Barcodes', y='Bases') + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
  data_long_pct = melt(data_pct, id.var = "READ_GROUP")
  my_sequences = levels(factor(data_long_pct$READ_GROUP))
  data_long_pct$READ_GROUP <- factor(data_long_pct$READ_GROUP, levels = data_long_pct$READ_GROUP)
  p2 = ggplot(data_long_pct, aes(x=READ_GROUP, y = value, fill = variable)) + geom_bar(stat = 'identity') + theme(axis.text.x=element_text(angle = 90, hjust = 0)) + labs(x='Barcodes', y='%Bases')
  # This allows to align the main plots so that we can relate both directly with the label from the bottom one.
  gp1 <- ggplotGrob(p1)
  gp2 <- ggplotGrob(p2)
  pdf(file = file.path(path,'plots',paste0(sample_name, '_rna_metrics.pdf')), width = 16, height = 13)
  grid::grid.newpage()
  grid::grid.draw(rbind(gp1, gp2))
  dev.off()
}


plotPolyTrim = function(file_path, sample_name, path){
  data = read.table(file = file_path, header=T, stringsAsFactors=F, skip = 4)
  polya = ggplot(data, aes(x=BIN, y = VALUE)) + geom_bar(stat = 'identity') + labs(title=paste('Length of poly trimmed in ', sample_name), x='Length', y='Counts') + theme_pubr()
  ggsave(plot = polya, paste0(sample_name, '_polya_trimmed.pdf'), path = paste0(path,'/plots/'))
}

for(i in 1:length(samples)){
  top_barcodes = config_file_data$Samples[[i]]$expected_cells
  plotRNAMetrics(file_path = paste0(path,"logs/",samples[i],"_rna_metrics.txt"), sample_name = samples[i], top_barcodes = top_barcodes, path = path)
  plotPolyTrim(file_path = paste0(path,"logs/",samples[i],"_polyA_trim.txt"), sample_name = samples[i], path = path)
}