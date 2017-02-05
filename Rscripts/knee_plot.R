library(jsonlite)
args = commandArgs(TRUE)

# You can tweak here the fraction of reads you expect
# when a cell is captured. I've found that this value often
# corresponds to a good minimum to validate a STAMP
fraction = 0.001

#You can tweak the multiplier of barcodes you need here
# to change the value on the X axis
xlim = 2.5 * fromJSON(paste0(path,'/config.json'))$Barcodes
path = args[1]

samples = names(fromJSON(paste0(path,'/config.json'))$Samples)

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
  pdf(file = paste0(path,"plots/",samples[i],"_knee_plot.pdf"), width = 5, height = 5)
  temp = plotCumulativePlot(file_path = paste0(path,samples[i],"_hist_out_cell.txt"), title = paste0(samples[i],"\nMinCellFraction = ", fraction), fraction = fraction, x_scale = xlim)
  dev.off()
  write.table(temp, file = paste0(path,"summary/",samples[i],"_barcodes.csv"),col.names = F, quote = F, row.names = F)
}
