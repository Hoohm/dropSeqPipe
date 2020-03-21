#------------------------------------ for debugging:
# For debugging add the following line in config.yaml (without the #)
# DEBUG: True
# This will create R objects in the debug directory containing the snakemake
# object R object that can be loaded into a custom R session as below:
debug_flag <- FALSE
if (snakemake@config$DEBUG) {
  debug_flag <- TRUE
  message("In debug mode: saving R objects to inspect later")
  path_debug <- file.path(snakemake@config$LOCAL$results, "debug")
  dir.create(path_debug)
  save(snakemake, file = file.path(path_debug, "plot_yield_snakemake.rdata"))
}
#------------------------------------ debugging

library(ggplot2)
library(tidyr)
library(grid)
library(gridExtra)
library(viridis)
library(stringr)

samples <- snakemake@params$sample_names
batches <- snakemake@params$batches
mydata  <- data.frame(matrix(nrow = length(samples), ncol = 7))
colnames(mydata) <- c("Sample", "Batch", "Cutadapt filtered", "Unmapped",
                      "Multi mapped", "Uniquely mapped", "Total reads")
mydata[, "Sample"] <- samples
mydata[, "Batch"]  <- batches
for (i in 1:length(samples)) {
  # Input files and variables
  STAR_output <- read.table(snakemake@input$STAR_output[i],
                            skip = 5, sep = "\t",
                            fill = TRUE, stringsAsFactors = FALSE)
  mysample <- samples[i]
  # Read files
  # bbmap_log = read.table(snakemake@input$repaired[i], sep=':', header=FALSE, skip=8, row.names=1, nrows=4)
  # reads_after_filtering = as.numeric(str_match(bbmap_log['Pairs',], pattern = "\t([0-9]{1,20}) reads.*")[,2])/2
  bbmap_log   <- readLines(snakemake@input$repaired[i])
  reads_after_filtering <- as.numeric(str_match(bbmap_log[grep("^Pairs:", bbmap_log)],
                                                pattern = "\t([0-9]{1,20}) reads.*")[, 2]) / 2
  R1_filtered <- read.table(snakemake@input$R1_filtered[i], header = FALSE, skip = 7, sep = ":", nrows = 7, row.names = 1)
  total_reads <- as.numeric(str_replace_all(R1_filtered["Total reads processed", ], pattern = (" |,"), ""))

  # R2_filtered = read.table(snakemake@input$R2_filtered[i], header = FALSE, skip=8, sep=':', nrows=7, row.names=1)

  # R1_adapters = as.numeric(str_remove_all(str_match(R1_filtered['Reads with adapters',], pattern = "(.*) \\(")[,2], pattern = (' |,')))
  # R1_too_short = as.numeric(str_remove_all(str_match(R1_filtered['Reads that were too short',], pattern = "(.*) \\(")[,2], pattern = (' |,')))
  # R1_passed = as.numeric(str_remove_all(str_match(R1_filtered['Reads written (passing filters)',], pattern = "(.*) \\(")[,2], pattern = (' |,')))
  # R1_filtered = total_reads - R1_passed

  # R2_adapters = as.numeric(str_remove_all(str_match(R2_filtered['Reads with adapters',], pattern = "(.*) \\(")[,2], pattern = (' |,')))
  # R2_too_short = as.numeric(str_remove_all(str_match(R2_filtered['Reads that were too short',], pattern = "(.*) \\(")[,2], pattern = (' |,')))
  # R2_passed = as.numeric(str_remove_all(str_match(R2_filtered['Reads written (passing filters)',], pattern = "(.*) \\(")[,2], pattern = (' |,')))
  # R2_filtered = total_reads - R2_passed

  mydata[which(mydata$Sample == mysample), "Cutadapt filtered"] <- total_reads - reads_after_filtering
  mydata[which(mydata$Sample == mysample), "Total reads"] <- total_reads

  # STAR output
  reads_in        <- as.numeric(STAR_output$V2[1])
  uniquely_mapped <- as.numeric(STAR_output$V2[4])
  multi_mapped    <- as.numeric(STAR_output$V2[19])
  unmapped        <- reads_in - uniquely_mapped - multi_mapped

  mydata[which(mydata$Sample == mysample), "Uniquely mapped"] <- uniquely_mapped
  mydata[which(mydata$Sample == mysample), "Multi mapped"]    <- multi_mapped
  mydata[which(mydata$Sample == mysample), "Unmapped"]        <- unmapped
}

# tidyr version
 mydata_long <- mydata %>% gather(variable, value, -Sample, -Batch)
# melt will be retired, use gather instead: https://github.com/hadley/reshape
#Force factor order.
mydata_long$variable = factor(mydata_long$variable, levels = c('Cutadapt filtered','Multi mapped','Total reads','Unmapped','Uniquely mapped'))
color_palette = c('#e88270','#cb7262','#ae6254','#70d6e8')


p1 <- ggplot(subset(mydata_long, mydata_long$variable != "Total reads"),
             aes(x = Sample, y = value, fill = variable)) +
  geom_histogram(stat = "identity", binwidth = 1 / length(samples)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  labs(title = paste("Yield of all the reads for each category"),
       x = "Samples",
       y = "Number of reads") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold")) +
  facet_grid(~Batch, scales = "free") +
  scale_fill_viridis(discrete = TRUE, option = "viridis") +
  scale_y_continuous(labels = scales::scientific)

mydata_pct <- mydata[, -c(1, 2)] / mydata[, "Total reads"]
mydata_pct <- cbind(Sample = mydata[, "Sample"],
                    Batch = mydata[, "Batch"], mydata_pct)

mydata_long_pct <- mydata_pct %>% gather(variable, value, -Sample, -Batch)

mydata_long_pct$variable = factor(mydata_long$variable, levels = c('Cutadapt filtered','Multi mapped','Unmapped','Uniquely mapped'))

p2 <- ggplot(subset(mydata_long_pct, mydata_long_pct$variable != "Total reads"),
             aes(x = Sample, y = value, fill = variable)) +
  labs(fill = "Filters") +
  geom_histogram(stat = "identity", binwidth = 1 / length(samples)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(x = "Samples",
       y = "Percentage of reads") +
  facet_grid(~Batch, scales = "free") +
  scale_fill_viridis(discrete = TRUE, option = "viridis") +
  scale_y_continuous(labels = scales::percent)

gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)

pdf(file = snakemake@output$pdf, width = 16, height = 13)
grid::grid.newpage()
grid::grid.draw(rbind(gp1, gp2, size = "last"))
dev.off()

if (debug_flag) {
  save.image(file = file.path(path_debug, "plot_yield_workspace.rdata"))
}
