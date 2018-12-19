library(ggplot2)
library(plyr)
# Create the cumulative plot
data=read.table(file = snakemake@input[[1]][1], header=FALSE, stringsAsFactors=FALSE)
barcodes = data$V2
total_reads = sum(data$V1)
y_raw=cumsum(data$V1)
y=(y_raw/total_reads)
x = 1:length(y)
plot_data = data.frame(rank = x,cum_sum=y, Barcode=data$V2)
x_scale = snakemake@params$cells * 4
knee_plot = ggplot(plot_data, aes(x=rank, y=cum_sum))
knee_plot = knee_plot + geom_point(size = 0.1) 
knee_plot = knee_plot + xlim(0,x_scale)
knee_plot = knee_plot + geom_vline(xintercept=snakemake@params$cells, linetype="dashed", color = "red")
knee_plot = knee_plot + ggtitle(paste0(snakemake@wildcards$sample, '\nTotal reads: ', prettyNum(total_reads)))
knee_plot = knee_plot + theme(plot.title = element_text(size=10))
knee_plot = knee_plot + labs(x='STAMPS', y='Cumulative fraction of reads')
knee_plot = knee_plot + scale_y_continuous(labels = scales::percent)

if(!is.null(snakemake@input$barcodes))
{
	selected_cells = read.csv(snakemake@input$barcodes, header=FALSE, stringsAsFactors=FALSE)
	knee_plot = knee_plot + geom_point(data = plot_data[plot_data$Barcode %in% selected_cells$V1,], aes(x=rank, y=cum_sum, color='Selected'), size=0.1)
	knee_plot = knee_plot + scale_color_manual(values=c('Selected'='green'))
 }
ggsave(knee_plot, file=snakemake@output$pdf, width = 4, height = 3)
