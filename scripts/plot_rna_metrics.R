library(ggplot2)
library(tidyr)
library(gridExtra)
library(grid)

data <- read.csv(file = snakemake@input$rna_metrics, header = T, stringsAsFactors = F, skip = 6, sep = "\t")
data <- data[order(data$PF_ALIGNED_BASES, decreasing = T), ]
data_pct <- data[, c("READ_GROUP", "PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES")]
data <- data[, c("READ_GROUP", "RIBOSOMAL_BASES", "CODING_BASES", "UTR_BASES", "INTRONIC_BASES", "INTERGENIC_BASES")]
data_long <- data %>% gather(read_overlap, count, -READ_GROUP)
# Keep the order of the barcodes using factor and levels.
my_sequences <- factor(unique(data_long$READ_GROUP))
data_long$READ_GROUP <- factor(data_long$READ_GROUP, levels = my_sequences)
p1 <- ggplot(data_long, aes(x = READ_GROUP, y = count, fill = read_overlap)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0))
if (!is.null(snakemake@params$cells)) {
  p1 <- p1 + labs(title = paste("Top", snakemake@params$cells, "barcodes for", snakemake@wildcards$sample),
                  x = "Barcodes", y = "Bases")
} else {
  p1 <- p1 + labs(title = paste("Whitelisted barcodes for", snakemake@wildcards$sample),
                  x = "Barcodes", y = "Bases")
}
p1 <- p1 + theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank())
p1 <- p1 + scale_y_continuous(labels = scales::scientific)

data_long_pct <- data_pct %>% gather(read_overlap, fraction, -READ_GROUP)
data_long_pct$READ_GROUP <- factor(data_long_pct$READ_GROUP,
                                   levels = my_sequences)
p2 <- ggplot(data_long_pct, aes(x = READ_GROUP, y = fraction, fill = read_overlap)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  labs(x = "Barcodes", y = "%Bases") +
  scale_y_continuous(labels = scales::percent)
# This allows to align the main plots so that we can relate both directly with the label from the bottom one.
gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)
pdf(file = snakemake@output$pdf, width = 16, height = 13)
grid::grid.newpage()
grid::grid.draw(rbind(gp1, gp2, size = "last"))
dev.off()
