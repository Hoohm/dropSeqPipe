library(yaml)
args = commandArgs(TRUE)
path= args[1]
samples = names(yaml.load_file(paste0(path,'config.yaml'))$Samples)

read_list = c()
for (i in 1:length(samples)){
  temp_matrix = read.table(paste0(path,"summary/",samples[i],'_expression_matrix.txt'), header=T, stringsAsFactors = F)
  old_names = colnames(temp_matrix)[-1]
  colnames(temp_matrix) = c("GENE",paste(samples[i], old_names, sep = "_"))
  read_list=c(read_list, list(temp_matrix))
}

# Little function that allows to merge unequal matrices
merge.all <- function(x, y) {
  merge(x, y, all=TRUE, by="GENE")
}

read_counts <- Reduce(merge.all, read_list)
read_counts[is.na(read_counts)] <- 0
rownames(read_counts) = read_counts[,1]
read_counts = read_counts[,-1]

write.table(read_counts, file=paste0(path,"summary/gene_name_matrix.tsv"), sep='\t')    