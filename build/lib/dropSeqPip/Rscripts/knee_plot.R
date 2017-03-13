library(yaml)
args = commandArgs(TRUE)
path = args[1]
config_file_data = yaml.load_file(paste0(path,'config.yaml'))

samples = names(config_file_data$Samples)

plotCumulativePlot = function(file_path, title, fraction, x_scale){
  data=read.table(file = file_path, header=F, stringsAsFactors=F)
  x_raw=cumsum(data$V1)
  x=x_raw/max(x_raw)
  plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",
       ylab="cumulative fraction of reads", xlim=c(1,x_scale), main=paste(title,"\nTotal reads: ",prettyNum(sum(as.numeric(x_raw)))))
  id = length(subset(diff(x)[seq(x_scale)],diff(x)[seq(x_scale)] > fraction))
  abline(v = id, col='red')
  abline(h = x[id], col='red')
  return(data$V2[seq(id)])
}


for(i in 1:length(samples)){
  pdf(file = paste0(path,"/plots/",samples[i],"_knee_plot.pdf"), width = 5, height = 5)
  fraction = config_file_data$Samples[[i]]$fraction
  temp = plotCumulativePlot(file_path = paste0(path,"summary/",samples[i],"_hist_out_cell.txt"), title = paste0(samples[i],"\nMinCellFraction = ", fraction), fraction = fraction, x_scale = config_file_data$Samples[[i]]$expected_cells * 5)
  dev.off()
  if(length(config_file_data$SPECIES) == 1){
    write.table(temp, file = paste0(path,"summary/", samples[i], '_', config_file_data$SPECIES, "_barcodes.csv"),col.names = F, quote = F, row.names = F)
  } else if (length(config_file_data$SPECIES) == 2){
    write.table(temp, file = paste0(path,"summary/",samples[i],"_barcodes.csv"),col.names = F, quote = F, row.names = F)
  }
}
