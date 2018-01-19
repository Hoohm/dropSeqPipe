library(ggplot2)
library(reshape2)
samples = snakemake@params$sample_names
data = data.frame(matrix(nrow=length(samples), ncol=511))
colnames(data) = c('Sample', seq(1:50))
data[,'Sample'] = samples

for(i in 1:length(samples)){
  sample = samples[i]
  umi_per_gene = read.csv(snakemake@input[[i]][1], sep='\t')
  freq_table = data.frame(table(umi_per_gene$Num_Obs))
  data[i,2:51] = as.numeric(freq_table$Freq[1:50])
}

data_long = melt(data, 'Sample')
data_long$variable = as.numeric(data_long$variable)
p = ggplot(data_long, aes(x=variable, y=value, color= Sample))
p = p + geom_line()
p = p + scale_y_log10()
p = p + scale_x_continuous(breaks = c(0,25,50))
p = p + ggtitle('UMI per gene distribution')
ggsave(plot = p, snakemake@output$pdf, height = 4, width = 6)
ggsave(plot = p, snakemake@output$png, height = 4, width = 6)