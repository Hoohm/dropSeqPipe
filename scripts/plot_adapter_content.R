library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(viridis)
library(stringr)





samples = snakemake@params$sample_names
batches = snakemake@params$batches
temp = read.csv(snakemake@input[[1]][1], header = TRUE)
data = data.frame(matrix(nrow=length(samples)*nrow(temp), ncol=5))
colnames(data) = c('Sample','Batch','Adapter','Pair','Count')

data[,'Sample'] = rep(samples, nrow(temp))
data[,'Batch'] = rep(batches,nrow(temp))



for(i in 1:length(samples)){
  sample = samples[i]
  #Read files
  cutadapt_clean = read.csv(snakemake@input[[i]][1], header = TRUE)
   
  data[which(data$Sample==sample),c(3,4,5)] = cutadapt_clean[,c('Adapter','Pair','Count')]
}

data$Adapter = factor(data$Adapter)
data$Pair = factor(data$Pair)
levels(data$Adapter) = levels(temp$Adapter)
levels(data$Pair) = levels(temp$Pair)

p1 = ggplot(data, aes(x=Sample, y = Count, fill = Adapter))
p1 = p1 + geom_bar(stat = 'identity')
p1 = p1 + theme(axis.text.x=element_text(angle = 90, hjust = 0))
p1 = p1 + facet_grid(~Batch, scales = "free")
p1 = p1 + facet_wrap(~Pair, nrow=2, scales="free") + theme_minimal()
p1 = p1 + ggtitle('Comparison accross samples of adapter content')
p1

ggsave(p1, snakemake@output$pdf)