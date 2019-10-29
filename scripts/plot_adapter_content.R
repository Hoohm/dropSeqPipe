#------------------------------------ for debugging:
# For debugging add the following line in config.yaml (without the #)
# DEBUG: True
# This will create R objects in the debug directory containing the snakemake
# object R object that can be loaded into a custom R session as below:

debug_flag <- FALSE
# check if DEBUG flag is set
if (snakemake@config$DEBUG) {
  debug_flag <- TRUE
  message("In debug mode: saving R objects to inspect later")
  path_debug <- file.path(snakemake@config$LOCAL$results, "debug")
  dir.create(path_debug, showWarnings = FALSE)
  save(snakemake, file = file.path(path_debug, "plot_adapter_content_snakemake.rdata"))
}

#------------------------------------ debugging

library(ggplot2)
library(dplyr)
library(viridis)

samples <- snakemake@params$sample_names
batches <- snakemake@params$batches

#Read files into a list
cutadapt_clean_list <- list()
for (i in seq_along(samples)){
  cutadapt_clean           <- read.csv(snakemake@input[[i]][1], header = TRUE)
  cutadapt_clean$Sample    <- samples[i]
  cutadapt_clean$Batch     <- batches[i]
  cutadapt_clean_list[[i]] <- cutadapt_clean
}

# combining adaptors accross samples
cutadapt_counts <- Reduce(rbind, cutadapt_clean_list, NULL)

#Transform it into percentages
cutadapt_counts <- group_by(cutadapt_counts, Sample, Pair) %>%
  mutate(Percentages=Count/sum(Count))
#      Adapter                           Sequence Pair Count  Sample  Batch
# 1 PrefixNX/1                AGATGTGTATAAGAGACAG   R1     7 sample1 Batch1
# ...
# 6  Trans2_rc CTGTCTCTTATACACATCTCCGAGCCCACGAGAC   R2     5 sample2 Batch2

p1 <- ggplot(cutadapt_counts, aes(x=Sample, y = Percentages, fill = Adapter))  +
  geom_bar(stat = "identity") +
  facet_grid(Pair ~ Batch, scales = "free") +
  theme_minimal() +
  ggtitle("Comparison accross samples of adapter content") +
  scale_x_discrete(label=abbreviate) +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  scale_fill_viridis(discrete=TRUE)

ggsave(plot=p1, filename=snakemake@output$pdf)

if (debug_flag) {
  save.image(file = file.path(path_debug, "plot_adapter_content_workspace.rdata"))
}
