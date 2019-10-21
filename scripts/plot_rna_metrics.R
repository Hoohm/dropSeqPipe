library(ggplot2)
library(tidyr)
library(gridExtra)
library(grid)
library(viridis)
debug_flag <- FALSE
if (snakemake@config$DEBUG) {
  debug_flag <- TRUE
  message("In debug mode: saving R objects to inspect later")
  path_debug <- file.path(snakemake@config$LOCAL$results, "debug")
  dir.create(path_debug, showWarnings = FALSE)
  save(snakemake, file = file.path(path_debug, "plot_rna_metrics_snakemake.rdata"))
}

#### /debug

mydata <- read.csv(file = snakemake@input$rna_metrics, header = T,
                   stringsAsFactors = F, skip = 6, sep = "\t")
mydata <- mydata[order(mydata$PF_ALIGNED_BASES, decreasing = T), ]
mydata_pct <- mydata[, c("READ_GROUP", "PCT_RIBOSOMAL_BASES",
                         "PCT_CODING_BASES", "PCT_UTR_BASES",
                         "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES")]
mydata <- mydata[, c("READ_GROUP", "RIBOSOMAL_BASES", "CODING_BASES",
                     "UTR_BASES", "INTRONIC_BASES", "INTERGENIC_BASES")]

# converting into long format for ploting
mydata_long <- mydata %>% gather(read_overlap, count, -READ_GROUP)

# Keep the original order of the barcodes using factor and levels.
mydata_long$READ_GROUP <- factor(mydata_long$READ_GROUP,
                                 levels = factor(unique(mydata_long$READ_GROUP)))
mydata_long$read_overlap <- factor(mydata_long$read_overlap,
                                   levels = unique(mydata_long$read_overlap))

p1 <- ggplot(mydata_long, aes(x = READ_GROUP, y = count, fill = read_overlap)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0))
p1 <- p1 + labs(title = paste(nrow(mydata),
                              "selected barcodes for",
                              snakemake@wildcards$sample),
                  x = "Barcodes", y = "Bases")
p1 <- p1 + theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank())
p1 <- p1 + scale_y_continuous(labels = scales::scientific)
p1 <- p1 + scale_fill_viridis(discrete = TRUE, option = "viridis")


mydata_long_pct <- mydata_pct %>% gather(read_overlap, fraction, -READ_GROUP)
# Keep the original order of the barcodes using factor and levels.
mydata_long_pct$READ_GROUP <- factor(mydata_long_pct$READ_GROUP,
                                 levels = factor(unique(mydata_long_pct$READ_GROUP)))
mydata_long_pct$read_overlap <- factor(mydata_long_pct$read_overlap,
                                   levels = unique(mydata_long_pct$read_overlap))

p2 <- ggplot(mydata_long_pct, aes(x = READ_GROUP, y = fraction, fill = read_overlap)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  labs(x = "Barcodes", y = "%Bases") +
  scale_y_continuous(labels = scales::percent) + scale_fill_viridis(discrete = TRUE, option = "viridis")
# This allows to align the main plots so that we can relate both directly with the label from the bottom one.
gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)
pdf(file = snakemake@output$pdf, width = 16, height = 13)
grid::grid.newpage()
grid::grid.draw(rbind(gp1, gp2, size = "last"))
dev.off()

if (debug_flag) {
  save.image(file = file.path(path_debug, "plot_rna_metrics_workspace.rdata"))
}
