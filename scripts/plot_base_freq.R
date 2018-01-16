library(ggplot2)
library(ggseqlogo)
# Create the frequency of nucleotide plot for cell barcodes

barcodes = read.table(file = snakemake@input[[1]], header = F, stringsAsFactors = F)
freq_plot = ggplot() + geom_logo(barcodes, method = 'probability', col_scheme = 'nucleotide')
freq_plot = freq_plot + ggtitle(snakemake@wildcards$sample)
ggsave(, file = snakemake@output[[1]], width = 4, height = 3)