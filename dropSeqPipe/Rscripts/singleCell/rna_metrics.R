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
  data_pct = data[,c("READ_GROUP","PCT_RIBOSOMAL_BASES","PCT_CODING_BASES","PCT_UTR_BASES","PCT_INTRONIC_BASES","PCT_INTERGENIC_BASES")]
  data = data[,c("READ_GROUP","RIBOSOMAL_BASES","CODING_BASES","UTR_BASES","INTRONIC_BASES","INTERGENIC_BASES")]
  data_long = melt(data, id.var = "READ_GROUP")
  # Keep the order of the barcodes using factor and levels.
  my_sequences = factor(unique(data_long$READ_GROUP))
  data_long$READ_GROUP <- factor(data_long$READ_GROUP, levels = my_sequences)
  p1 = ggplot(data_long, aes(x=READ_GROUP, y = value, fill = variable)) + geom_bar(stat = 'identity') + theme(axis.text.x=element_text(angle = 90, hjust = 0)) + labs(title=paste('Top',top_barcodes, 'barcodes for', sample_name), x='Barcodes', y='Bases') + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
  data_long_pct = melt(data_pct, id.var = "READ_GROUP")
  data_long_pct$READ_GROUP <- factor(data_long_pct$READ_GROUP, levels = my_sequences)
  p2 = ggplot(data_long_pct, aes(x=READ_GROUP, y = value, fill = variable)) + geom_bar(stat = 'identity') + theme(axis.text.x=element_text(angle = 90, hjust = 0)) + labs(x='Barcodes', y='%Bases') + scale_y_continuous(labels = scales::percent)
  # This allows to align the main plots so that we can relate both directly with the label from the bottom one.
  gp1 <- ggplotGrob(p1)
  gp2 <- ggplotGrob(p2)
  pdf(file = file.path(path,'plots',paste0(sample_name, '_rna_metrics.pdf')), width = 16, height = 13)
  grid::grid.newpage()
  grid::grid.draw(rbind(gp1, gp2))
  dev.off()
}

plotBCDrop = function(path, samples, config_file_data){
  data = data.frame(matrix(nrow=length(samples), ncol=4))
  data[,4] = samples
  colnames(data) = c('BC_drop','UMI_drop','Total_reads','Sample')
  fastqc_summary = read.table(file.path(path, 'summary', 'fastqc.txt'), header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  BC_length = config_file_data$GLOBAL$Cell_barcode$end-config_file_data$GLOBAL$Cell_barcode$start +1
  UMI_length = config_file_data$GLOBAL$UMI$end-config_file_data$GLOBAL$UMI$start +1
  for(i in 1:length(samples)){
    sample = samples[i]
    data[i,3] = as.integer(unique(as.character(fastqc_summary['Total Sequences',grep(sample, colnames(fastqc_summary))])))
    
    BC_drop_data = read.table(file.path(path, 'logs',paste0(sample,'_CELL_barcode.txt')), header = TRUE, row.names = 1)
    num_reads_BC_dropped = sum(BC_drop_data[(config_file_data$GLOBAL$Cell_barcode$num_below_quality+1):(1+BC_length),])
    data[i,1] = num_reads_BC_dropped
    
    UMI_drop_data = read.table(file.path(path, 'logs',paste0(sample,'_UMI_barcode.txt')), header = TRUE, row.names = 1)
    num_reads_UMI_dropped = sum(UMI_drop_data[(config_file_data$GLOBAL$UMI$num_below_quality+1):(1+UMI_length),])
    data[i,2] = num_reads_UMI_dropped
  }
  
  data_long = melt(data, 'Sample')
  # Keep the order of the barcodes using factor and levels.
  p1 = ggplot(data_long, aes(x=Sample, y = value, fill = variable)) + geom_bar(stat = 'identity') + theme(axis.text.x=element_text(angle = 90, hjust = 0)) + labs(title=paste('BC and UMI drop comparison.\nMin Cell BC quality =',config_file_data$GLOBAL$Cell_barcode$min_quality,'& Max number bellow =',config_file_data$GLOBAL$Cell_barcode$num_below_quality,'\nMin UMI BC quality =',config_file_data$GLOBAL$UMI$min_quality,'& Max number bellow =',config_file_data$GLOBAL$UMI$num_below_quality), x='Samples', y='Number of reads') + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())  +scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  data$BC_drop_pct = data$BC_drop/data$Total_reads
  data$UMI_drop_pct = data$UMI_drop/data$Total_reads
  data$not_dropped = (data$Total_reads-data$BC_drop-data$UMI_drop)/data$Total_reads
  data_long_pct = melt(data,'Sample', measure.vars = c('BC_drop_pct', 'UMI_drop_pct'))
  p2 = ggplot(data_long_pct, aes(x=Sample, y = value, fill = variable)) + geom_bar(stat = 'identity') + theme(axis.text.x=element_text(angle = 90, hjust = 0)) + labs(x='Samples', y='Percentage of reads')+ scale_y_continuous(labels = scales::percent) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  # This allows to align the main plots so that we can relate both directly with the label from the bottom one.
  gp1 <- ggplotGrob(p1)
  gp2 <- ggplotGrob(p2)
  pdf(file = file.path(path,'plots','BCDrop.pdf'), width = 8, height = 7)
  grid::grid.newpage()
  grid::grid.draw(rbind(gp1, gp2))
  dev.off()
}


plotPolyATrim = function(file_path, sample_name, path){
  data = read.table(file = file_path, header=T, stringsAsFactors=F, skip = 4)
  polya = ggplot(data, aes(x=BIN, y = VALUE)) + geom_bar(stat = 'identity') + labs(title=paste('Length of polyA trimmed in\n', sample_name), x='Length', y='Counts') + theme_pubr() +  geom_smooth() + scale_x_continuous(breaks = seq(0,data$BIN[length(data$BIN)], 10), labels = seq(6,data$BIN[length(data$BIN)]+6, 10))
  ggsave(plot = polya, paste0(sample_name, '_polya_trimmed.pdf'), path = paste0(path,'/plots/'), height = 6, width = 8)
}

plotStartTrim = function(file_path, sample_name, path){
  data = read.table(file = file_path, header=T, stringsAsFactors=F, skip = 4)
  polya = ggplot(data, aes(x=BIN, y = VALUE)) + geom_bar(stat = 'identity') + labs(title=paste('Length of SMART adapter trimmed in\n', sample_name), x='Length', y='Counts') + theme_pubr() + scale_x_continuous(breaks = c(data$BIN), labels = factor(data$BIN))
  ggsave(plot = polya, paste0(sample_name, '_start_trimmed.pdf'), path = paste0(path,'/plots/'), height = 6, width = 8)
}

plotBCQuality = function(file_path, sample_name, path, config_file_data){
  data = read.table(file = file_path, header=T, stringsAsFactors=F)
  polya = ggplot(data[-1,], aes(x=num_failed_bases, y=num_barcodes)) + geom_bar(stat = 'identity') + labs(title=paste('Number of Cell barcode bases under',config_file_data$GLOBAL$Cell_barcode$min_quality,'quality in\n', sample_name), x='Num of failed bases', y='Counts') + theme_pubr()
  ggsave(plot = polya, paste0(sample_name, '_BC_qual.pdf'), path = paste0(path,'/plots/'), height = 4, width = 6)
}

plotUMIQuality = function(file_path, sample_name, path, config_file_data){
  data = read.table(file = file_path, header=T, stringsAsFactors=F)
  polya = ggplot(data[-1,], aes(x=num_failed_bases, y=num_barcodes)) + geom_bar(stat = 'identity') + labs(title=paste('Number of UMI barcode bases under', config_file_data$GLOBAL$UMI$min_quality,'quality in\n', sample_name), x='Num of failed bases', y='Counts') + theme_pubr()
  ggsave(plot = polya, paste0(sample_name, '_UMI_qual.pdf'), path = paste0(path,'/plots/'), height = 4, width = 6)
}

for(i in 1:length(samples)){
  top_barcodes = config_file_data$Samples[[i]]$expected_cells
  print('RNAmetrics')
  plotRNAMetrics(file_path = paste0(path,"logs/",samples[i],"_rna_metrics.txt"), sample_name = samples[i], top_barcodes = top_barcodes, path = path)
  plotPolyATrim(file_path = paste0(path,"logs/",samples[i],"_polyA_trim.txt"), sample_name = samples[i], path = path)
  plotStartTrim(file_path = paste0(path,"logs/",samples[i],"_start_trim.txt"), sample_name = samples[i], path = path)
  plotBCQuality(file_path = paste0(path,"logs/",samples[i],"_CELL_barcode.txt"), sample_name = samples[i], path = path,  config_file_data)
  plotUMIQuality(file_path = paste0(path,"logs/",samples[i],"_UMI_barcode.txt"), sample_name = samples[i], path = path, config_file_data)
}
plotBCDrop(path, samples, config_file_data)
