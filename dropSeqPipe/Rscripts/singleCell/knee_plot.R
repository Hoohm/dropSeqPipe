library(yaml)
library(ggplot2)
library(ggseqlogo)
args = commandArgs(TRUE)
path = args[1]
config_file_data = yaml.load_file(paste0(path,'config.yaml'))

samples = names(config_file_data$Samples)

# Create the cumulative plot
plotCumulativePlot = function(file_path, title, fraction, x_scale, path, sample_name){
  data=read.table(file = file_path, header=F, stringsAsFactors=F)
  barcodes = data$V2
  total_reads = sum(data$V1)
  y_raw=cumsum(data$V1)
  y=(y_raw/total_reads)
  x = 1:length(y)
  data = data.frame(cbind(x,y))
  id = length(subset(diff(y)[seq(x_scale)],diff(y)[seq(x_scale)] > fraction))
  knee_plot = ggplot(data, aes(x=x, y=y)) + geom_point(size = 0.1) +xlim(0,x_scale) + theme_minimal() + geom_vline(xintercept=id, linetype="dashed", color = "red") + geom_hline(yintercept=y[id], linetype="dashed", color = "red")+ ggtitle(paste0(title, '\nTotal reads: ', prettyNum(total_reads))) + theme(plot.title = element_text(size=10)) + labs(x='STAMPS', y='Cumulative fraction of reads')
  ggsave(plot=knee_plot, paste0(sample_name, '_knee_plot.pdf'), path = paste0(path,'/plots/'), width = 4, height = 3)
  return(barcodes[seq(id)])
}

# Create the frequency of nucleotide plot for cell barcodes
plotBaseFreqBarcode = function(file_path, sample_name, path){
  barcodes = read.table(file = file_path, header = F, stringsAsFactors = F)
  p2 = ggplot() + geom_logo(barcodes, method = 'probability', col_scheme = 'nucleotide') + ggtitle(sample_name)
  ggsave(plot = p2, paste0(sample_name, '_base_freq_plot.pdf'), path = paste0(path,'/plots/'), width = 4, height = 3)
}

for(i in 1:length(samples)){
  fraction = config_file_data$Samples[[i]]$fraction
  temp = plotCumulativePlot(file_path = paste0(path,"logs/",samples[i],"_hist_out_cell.txt"), title = paste0(samples[i],"\nMinCellFraction = ", fraction), fraction = fraction, x_scale = config_file_data$Samples[[i]]$expected_cells * 4, path = path, sample_name = samples[i])
  write.table(temp, file = paste0(path,"summary/",samples[i],"_barcodes.csv"),col.names = F, quote = F, row.names = F)
  plotBaseFreqBarcode(file_path = paste0(path,"summary/",samples[i],"_barcodes.csv"), sample_name = samples[i], path = path)
  
}
