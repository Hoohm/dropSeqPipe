#Functions used to plot the species plot for drop-seq mixed protocol
# Credit to James Nemesh

library(yaml)
args = commandArgs(TRUE)
path = args[1]

config_file_data = yaml.load_file(paste0(path,'/config.yaml'))
samples = config_file_data$Samples
species = config_file_data$SPECIES
#num_cells = fromJSON(paste0(path,'/config.yaml'))$Barcodes



categorizeCellsUsingKneeKnownNumCellsPaper<-function (digitalExpressionFileO1, digitalExpressionFileO2, organismOne, organismTwo, pureRatio=0.2, numCells, numBeads, point.cex=1.5, category='transcripts', xlim_range=NULL,ylim_range=NULL) {
  dfFull=getNumTranscriptsPerCellBarcodeByOrganismPair(digitalExpressionFileO1, digitalExpressionFileO2, organismOne, organismTwo, category)
  dfFull=dfFull[order(dfFull$total, decreasing=T),]
  dfFull$ratio_one=dfFull[,2]/dfFull[,4]
  dfFull=head (dfFull, n=numBeads)
  df=head (dfFull, n=numCells)
  
  dfNoCall=dfFull[-1:-numCells,]
  if (dim(dfNoCall)[1]>0) {
    dfNoCall$organism="No Call"   
  }
  
  df$organism="Mixed"
  
  idx=which(df$ratio_one>= (1-pureRatio))
  df[idx,]$organism=organismOne
  idx=which(df$ratio_one<= (pureRatio))
  df[idx,]$organism=organismTwo
  
  result=rbind(df, dfNoCall)
  
  maxRange=max(result[,2], result[,3])
  
  dforganismOne= result[result$organism==organismOne,]
  dforganismTwo= result[result$organism==organismTwo,]
  dfMixed= result[result$organism=="Mixed",]
  dfNoCall = result[result$organism=="No Call",]
  
  if (is.null(xlim_range)) {
    xlim_range=c(0,maxRange)
  }
  
  if (is.null(ylim_range)) {
    ylim_range=c(0,maxRange)
  }
  colors=c("blue","red", "purple", "grey")
  plot(dforganismOne[,2], dforganismOne[,3], col= colors[1], pch=16, xlim=xlim_range, ylim=ylim_range, xlab=paste(organismOne,category), ylab=paste(organismTwo,category), cex=point.cex)
  points(dforganismTwo[,2], dforganismTwo[,3], col= colors[2], pch=16, cex=point.cex)
  points(dfMixed[,2], dfMixed[,3], col= colors[3], pch=16, cex=point.cex)
  points(dfNoCall[,2], dfNoCall[,3], col= colors[4], pch=16, cex=point.cex)
  l=c(paste(organismOne, dim(dforganismOne)[1]), paste(organismTwo, dim(dforganismTwo)[1]), paste("Mixed", dim(dfMixed)[1]), paste("No Call", dim(dfNoCall)[1]))
  legend("topright", legend=l, fill= colors)
  title(paste('Species plot based on',category))
  return (df)   
  
}

getNumTranscriptsPerCellBarcodeByOrganismPair<-function (digitalExpressionFileO1, digitalExpressionFileO2, organismOne, organismTwo, category) {
  if (is.null(organismOne) || is.null(organismTwo)) return(NULL)
  
  o1=getGenesAndTranscriptsPerCellBarcode(digitalExpressionFileO1)
  o2=getGenesAndTranscriptsPerCellBarcode(digitalExpressionFileO2)
  
  commonBC=union(o1$cellBC, o2$cellBC)
  o1p=o1[match(commonBC, o1$cellBC),]
  o2p=o2[match(commonBC, o2$cellBC),]
  if(category == 'genes'){
    df=data.frame(tag=commonBC, o1Count=o1p$numGenes, o2Count=o2p$numGenes, stringsAsFactors=F)
  }
  else {
    df=data.frame(tag=commonBC, o1Count=o1p$numTranscripts, o2Count=o2p$numTranscripts, stringsAsFactors=F)
  }
  
  idx1=which(is.na(df$o1Count))
  idx2=which(is.na(df$o2Count))
  if (length(idx1)>0) df[idx1,]$o1Count=0
  if (length(idx2)>0) df[idx2,]$o2Count=0
  
  df$total=apply(df[,2:3], 1, sum, na.rm=T)
  df=df[order(df$total, decreasing=T),]
  colnames(df)[2]= organismOne
  colnames(df)[3]= organismTwo
  return (df)
}


getGenesAndTranscriptsPerCellBarcode<-function (digitalExpressionFile) {
  a=read.table(digitalExpressionFile, header=T, stringsAsFactors=F)
  colnames(a)=c("cellBC", "numGenes", "numTranscripts")
  return (a)  
}

for (i in 1:length(samples)){
  digitalExpressionFileO1 = paste0(path, 'summary/', names(samples)[i],'_',species[[1]],'_dge.summary.txt')
  digitalExpressionFileO2 = paste0(path, 'summary/', names(samples)[i],'_',species[[2]],'_dge.summary.txt')
  num_cells = config_file_data$Samples[[i]]$expected_cells
  par(mar=c(5,4,4,2)+0.5)
  pdf(paste0(path,'plots/',names(samples)[i],'_species_plot_genes.pdf'), height=8, width=8)
  categorizeCellsUsingKneeKnownNumCellsPaper(digitalExpressionFileO1, digitalExpressionFileO2, species[[1]], species[[2]], pureRatio= 0.3, numCells= num_cells, numBeads = num_cells * 2, point.cex= 1, category = 'genes')
  dev.off()
  pdf(paste0(path,'plots/',names(samples)[i],'_species_plot_transcripts.pdf'), height=8, width=8)
  df=categorizeCellsUsingKneeKnownNumCellsPaper(digitalExpressionFileO1, digitalExpressionFileO2, species[[1]], species[[2]], pureRatio= 0.3, numCells= num_cells, numBeads = num_cells * 2, point.cex= 1, category = 'transcripts')
  dev.off()
  organism1 = subset(df, df$organism ==species[[1]])
  organism2 = subset(df, df$organism ==species[[2]])
  out_file_organism1 = paste0(path, 'summary/',names(samples)[i], '_',species[[1]],'_barcodes.csv')
  out_file_organism2 = paste0(path, 'summary/',names(samples)[i], '_',species[[2]],'_barcodes.csv')
  write.table(organism1$tag, out_file_organism1, row.names=F, col.names=F, quote=F)
  write.table(organism2$tag, out_file_organism2, row.names=F, col.names=F, quote=F)
}