#I am not the author of all of that code. I can't find it on github anymore. If you know who wrote it, please tell me.

library(reshape)
library(stringr)
library(grid)
library(plyr)
library(ggplot2)

args = commandArgs(TRUE)
folder_path = args[1]
qc.folders = file.path(folder_path, 'fastqc')
out.file = file.path(folder_path,'summary/fastqc')
fastqc_plot = file.path(folder_path,'plots/fastqc.pdf')

#===============================================================================
# FastQC DATA READ FUNCTIONS

# => Internal Function to collect Base statistics
read.base <- function(linn){
  base.stat <- matrix(NA, nrow =12 , ncol =2)
  np <- 0
  for(i in 1:length(linn)){
    if(grepl('^>>', linn[i]) & !grepl('^>>END_MODULE', linn[i])){
      np <- np+1
      sp.line <- unlist(strsplit(linn[i], split = '\t' ))
      base.stat[np,1] <- gsub('>>','',sp.line[1])
      base.stat[np,2] <- sp.line[2]
    }
  }
  base.stat <- data.frame(stats = base.stat[,1],
                          summary = base.stat[,2],
                          stringsAsFactors = FALSE)
  return(base.stat)
}

# Internal Function to read out data.frames
read.out <- function(linn, df.string, col.name){
  out <- list()
  rec = FALSE
  for(i in 1:length(linn)){
    if(grepl(df.string, linn[i])){
      rec <- TRUE
      rc <- -2
    }       
    if(grepl("^>>END_MODULE", linn[i]))
      rec <- FALSE
    if(rec){
      rc <- rc+1
      if(rc > 0){
        spl <- unlist( strsplit( linn[i], split ='\t' ))[1:2]
        out[[rc]] <- spl
      }
    }
  }
  out <- do.call(rbind, out)
  out <- data.frame(base = out[,1],
                    col.name = out[,2],
                    stringsAsFactors = TRUE)
  return(out)
}

# FUNCTION To get the data columns right
cleanup <- function(x){
  out <- ldply(x)
  ord <- as.numeric(sapply(strsplit(levels(as.factor(out$base)), '-'), '[[', 1))
  out$base <- factor(out$base, levels = levels(out$base)[order(ord)]  )
  out$col.name <- as.numeric(as.character(out$col.name))
  colnames(out) <- c('file','x','y')
  return( out )
}

# FUNCTION to seperate R1 and R2
sep.reads <- function(x, phrase){
  out <- list()
  for(i in names(x)){
    out[[i]] <- x[[i]][ grepl(phrase, names(x[[i]])) ]
  }
  return(out)
}

#===============================================================================
# PLOT FUNCTIONS



# plot theme
my.theme = theme_bw()+
  theme(axis.text.x = element_text(angle= 45, hjust=1, size = 5),
        #legend.position ="none",
        legend.text = element_text(size = 6 ),
        strip.text = element_text(vjust=0.5, hjust=0.5),
        plot.title = element_text(size=20, hjust=0.5))

# 1) Tile plot for the BASE.STATS
gg.base <- function(base.stats, title){
  title = paste(title, "SUMMARY")
  f.ord <- base.stats[[1]]$stats
  base.stats <- ldply(base.stats)
  base.stats$stats <- factor(base.stats$stats, levels = rev(f.ord))
  base.stats$summary <- factor(base.stats$summary, levels = c('pass','warn','fail'))
  #base.stats$.id <- gsub(' ', '-', 
  #                       str_wrap(gsub('_|-',' ',base.stats$.id), width = 15))
  # plot
  gg <- ggplot(base.stats)+aes(x = .id, y = stats, fill = summary)+
    geom_tile(colour = "white", lwd = 1)+
    labs(x="", y="", title = title)+
    scale_fill_manual(values = c("darkolivegreen3","darkorange", "firebrick2"))+
    my.theme
  return(gg)
}

# 2) Line plot for the BASE.QUAL
gg.qual <- function(seq.qual, title){
  title = paste(title, "Mean Base Quality")
  in.df2 <- cleanup(seq.qual)
  in.df2$file <- gsub(' ', '-', 
                      str_wrap(gsub('_|-',' ', in.df2$file), width = 25))
  # set up color background
  background <- data.frame(lower = c(0, 20, 28), 
                           upper = c(20, 28, 40),
                           col = letters[1:3],
                           xmax = length(levels(in.df2$x)))
  gg <- ggplot()+
    ylim(c(0,40))+
    geom_line(data = in.df2,
              aes(x = x , y = y, colour = file, group = file))+
    geom_rect(data = background,
              mapping = aes(xmin = 0, xmax = xmax, ymin = lower, ymax = upper, fill = col ),
              alpha = .2 , show.legend = FALSE)+
    geom_line(data = in.df2,
              aes(x = x, y = y, colour = file, group = file))+
    labs(x="", y="Mean Base Quality", title = title)+
    my.theme
  return(gg)
}

#===============================================================================

qcs <- list.dirs(path = qc.folders, full.names = TRUE, recursive = F)
lv <- vector('list', length(qcs))
names(lv) <- gsub("_fastqc",'', basename(qcs))

full.list <- list('base' =lv, 'qual' = lv, 'len' = lv, summ = lv )


for(i in 1:length(qcs)){
  # open file conn
  file.conn <- file.path(qcs[i], "fastqc_data.txt")
  conn <- file( file.conn, open = 'r' )
  linn <-readLines( conn )
  # read the base stats
  full.list[['summ']][[i]] <- read.out(linn, '>>Basic Statistics', 'Value')
  full.list[['base']][[i]] <- read.base(linn)
  full.list[['qual']][[i]] <- read.out(linn, '>>Per base sequence quality', 'mean.qual')
  # close conn
  close(conn)
}

# Separarte R1 and R2
r1 <- sep.reads(full.list, '_R1')
r2 <- sep.reads(full.list, '_R2')

# Make the grid plot
#--------------------
gg.linkF <- function(gfun, gdat, title){
  f1 <- do.call( paste0('gg.',gfun), list(gdat[[gfun]], title))
  f2 <- try( ggplot_gtable(ggplot_build( f1 ))) 
  return(f2)
}


# set the plotting size according to the number of samples
pdf(fastqc_plot,
    width = (13+ceiling(length(r1)/5)),
    height = 20)
# setup the grid

grid.newpage()

my.layout <- grid.layout(5,5, width = c(0.05,1,0.0,1,0.05),
                         height=c(0.1,0.7,1,1,0.1))
pushViewport(viewport(layout = my.layout))
gp = gpar(col="grey90", lwd=5)


# plot READ1 in Col 2
pushViewport( viewport(layout.pos.col=2, layout.pos.row= 2))
gg.link.plt <- gg.linkF( 'base', r1, 'READ1' )
grid.draw(gg.link.plt)
grid.roundrect(gp =gp)
upViewport()

pushViewport( viewport(layout.pos.col=2, layout.pos.row= 3:4))
gg.link.plt <- gg.linkF( 'qual', r1, 'READ1' )
grid.draw(gg.link.plt)
grid.roundrect(gp =gp)
upViewport()


# plot READ1 in Col 2
pushViewport( viewport(layout.pos.col=4, layout.pos.row= 2))
gg.link.plt <- gg.linkF( 'base', r2, 'READ2' )
grid.draw(gg.link.plt)
grid.roundrect(gp =gp)
upViewport()

pushViewport( viewport(layout.pos.col=4, layout.pos.row= 3:4))
gg.link.plt <- gg.linkF( 'qual', r2, 'READ2' )
grid.draw(gg.link.plt)
grid.roundrect(gp =gp)
upViewport()
dev.off()


# Make the summary table
f.txt <- ldply(full.list$summ)
f.txt <- cast(f.txt, base ~ .id, value = "col.name")
write.table(f.txt, file = paste(out.file,".txt", sep =""), sep = "\t", row.names = FALSE)