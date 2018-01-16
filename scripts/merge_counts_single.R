#print(snakemake)
samples = snakemake@params$sample_names
read_list = c()
for (i in 1:length(samples)){
	temp_matrix = read.table(snakemake@input[[i]][1], header=T, stringsAsFactors = F)
	cell_barcodes = colnames(temp_matrix)[-1]
	colnames(temp_matrix) = c("GENE",paste(samples[i], cell_barcodes, sep = "_"))
	read_list=c(read_list, list(temp_matrix))
}

# Little function that allows to merge unequal matrices
merge.all <- function(x, y) {
  merge(x, y, all=TRUE, by="GENE")
}

read_counts <- Reduce(merge.all, read_list)
read_counts[is.na(read_counts)] = 0
rownames(read_counts) = read_counts[,1]
read_counts = read_counts[,-1]
write.table(read_counts, file=snakemake@output[[1]], sep='\t')