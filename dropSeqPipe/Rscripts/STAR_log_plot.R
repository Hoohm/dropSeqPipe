#I am not the author of all of that code. I can't find it on github anymore. If you know who wrote it, please tell me.

library(ggplot2)
library(reshape)
library(gtable)
library(grid)
args = commandArgs(TRUE)
folder_path = args[1]


######################################################################
# process STAR Logs
######################################################################

reader <- function(df){
  data = read.table(df, sep ="\t", row.names=1,
                    skip =4, fill =TRUE)
  out = as.numeric(gsub("%", "", data[,1]))
  return(out)
}

trim.leading <- function (x)  sub("^\\s+", "", x)

star.logs <- function(path=NULL){
  if(is.null(path)){
    path=getwd()
  }
  logs = list.files(path=path, pattern= "final.out")
  logsF = file.path(path, logs)
  data = do.call("cbind", lapply(logsF, function(df) reader(df)))
  data = data.frame(data)
  colnames(data) = gsub("Log.final.out", "", logs)
  rd = read.table(logsF[1], sep ="\t", 
                  skip =4, fill =TRUE)[,1]
  rd = as.character(gsub(" \\|", "", rd ))
  rd = as.character(trim.leading(rd))
  row.names(data) = gsub("  ", "", rd)
  return(data)
}

# GGPLOT

# subplotting function
plot.gg1 <- function(x, main){
  my.theme = theme_bw()+
    theme(axis.text.x = element_text(angle= 90, vjust=0.5, hjust=1),
          legend.position ="none",
          strip.text = element_text(vjust=0.5, hjust=0.5),
          plot.title = element_text(size=20, hjust=0.5))
  #plot.margin= unit(c(0.5, 0.5, -0.5, 0.5), "lines")
  gg = ggplot(x)+aes(x=variable, y =value)+
    geom_bar(stat="identity",aes(fill =variable))+
    facet_grid(ID ~., scales="free_y")+
    xlab("")+ylab("")+ggtitle(main)+
    my.theme
  return(gg)
}

plot.gg2 <- function(x, main){
  my.theme = theme_bw()+
    theme(axis.text.x = element_text(angle= 90, vjust=0.5, hjust=1),
          legend.position ="none",
          strip.text = element_text(hjust = 0.5, vjust = 0.5),
          plot.title = element_text(size=20, hjust=0.5))
  #plot.margin= unit(c(0.5, 0.5, -0.5, 0.5), "lines")
  gg = ggplot(x)+aes(x=variable, y =value)+
    geom_bar(stat="identity",aes(fill =variable))+
    facet_wrap(~ID, scales="free_y")+
    xlab("")+ylab("")+ggtitle(main)+
    my.theme
  return(gg)
}

# main plotting function
plot.all <- function(x){
  d1 = melt(x[c(1,2),])
  d2 = melt(x[c(4,5,6,7),])
  d3 = melt(x[c(19,20),])
  d4 = melt(x[c(25,26),])
  # get gtables
  gt1 = ggplot_gtable(ggplot_build(plot.gg2(d1, main="Input reads")))
  gt2 = ggplot_gtable(ggplot_build(plot.gg1(d2, main="Uniquely mapped reads")))
  gt3 = ggplot_gtable(ggplot_build(plot.gg1(d3, main="Multi-mapping reads")))
  gt4 = ggplot_gtable(ggplot_build(plot.gg1(d4, main="Unmapped reads")))
  # gp
  gp = gpar(col="grey90", lwd=10)
  # gtable
  my.layout <- grid.layout(7,5, width = c(0.1,1,1,1,0.1),
                           height=c(0.1,1,1,1,1,1,0.1))
  pushViewport(viewport(layout = my.layout))
  pushViewport( viewport(layout.pos.col=(2), layout.pos.row=2:6))
  grid.draw(gt2)
  grid.roundrect(gp =gp)
  upViewport()
  pushViewport( viewport(layout.pos.col=(3:4), layout.pos.row=2:3))
  grid.draw(gt1)
  grid.roundrect(gp = gp)
  upViewport()
  pushViewport( viewport(layout.pos.col=(3), layout.pos.row=4:6))
  grid.draw(gt3)
  grid.roundrect(gp =gp)
  upViewport()
  pushViewport( viewport(layout.pos.col=(4), layout.pos.row=4:6))
  grid.draw(gt4)
  grid.roundrect(gp = gp)
  upViewport()
}

data <- star.logs(file.path(folder_path, 'logs'))
write.table(data, file = file.path(folder_path, "summary/STAR_Log_Stat.txt"), sep="\t", col.names=T, row.names=T, quote=F)
data$ID <- factor(row.names(data), levels = row.names(data))
pdf(file.path(folder_path, "plots/STAR_Log_Stat.pdf"), height=12, width=13)
plot.all(data)
dev.off()
