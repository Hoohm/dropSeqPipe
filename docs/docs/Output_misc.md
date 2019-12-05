# Other output

Apart from plots, dSP also generates various other output, such as logs, tables, R objects, reports and intermediated files.

## Tables
Various summary statistics are also computed such as `barcode_stats_pre_filter.csv` and `barcode_stats_post_filter.csv`:

###' barcode_stats_pre_filter.csv
Various statistic (Total_raw_reads, Nr_barcodes_total	Nr_barcodes_more_than_1_reads	Nr_barcodes_more_than_10_reads	percentile99	percentile95	percentile50	Gini-index	Reads_assigned_to_expected_STAMPs) 

### barcode_stats_post_filter.csv

Statistics for barcodes that were taken forward as STAMPs as set as `expected_cells` in `config.yaml`, including:
Total_nb_reads	Nb_STAMPS	Median_reads_per_STAMP	Mean_reads_per_STAMP	Total_nb_UMIs	Median_UMIs_per_STAMP	Mean_UMIs_per_STAMP	Mean_UMIs_per_Gene	Median_number_genes_per_STAMP	Mean_number_genes_per_STAMP	Mean_Ribo_pct	Mean_Mito_pct	Mean_Count_per_UMI	Read_length	Number_barcodes_used_for_debug	Pct_reads_after_filter_expected_cells	Pct_reads_after_filter_everything

## R objects

### Seurat 

A Seurat object is generated automatically which can be used straight away using Seurat and can be found here:
`summary/R_Seurat_objects.rdata`
Note, dSP up to version 0.5 produces Seurat v2 object. Newer version (dSP 0.6+) generate Seurat v3 objects.

As from dSP version 0.6 Seurat, you can import the object as follows:

'''
library(Seurat)
dspData <- readRDS("path/to/working_dir/results/summary/R_Seurat_objects.rdata")
'''
