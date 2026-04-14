#!/usr/bin/env Rscript

library(argparse)

# Parse command line arguments
parser <- ArgumentParser(description="OmniBenchmark module")

# Required by OmniBenchmark
parser$add_argument("--output_dir", dest="output_dir", type="character", required=TRUE,
                   help="Output directory for results")
parser$add_argument("--name", dest="name", type="character", required=TRUE,
                   help="Module name/identifier")
# Add your custom input arguments here
# Example:
# parser$add_argument("--input", dest="input", type="character", help="Input file")

# parameter for 'type' of filtering = ['manual', 'scrapper-auto']
parser$add_argument("--filter_type", dest="filter_type", type="character", help="type of filtering")

# parameter for 'input_h5' (comes from 1st stage)
parser$add_argument("--rawdata.h5ad", dest="input_h5", type="character", help="input file")

args <- parser$parse_args()

cat("Output directory:", args$output_dir, "\n")
cat("Module name:", args$name, "\n")

cat("Input file:", args$input_h5, "\n")
cat("Filtering type:", args$filter_type, "\n")

library(SingleCellExperiment)
library(anndataR)
library(scrapper)

sce <- read_h5ad(args$input_h5, as = "SingleCellExperiment")

is.mito <- grepl("^MT-", rownames(sce))
rna.qc.metrics <- computeRnaQcMetrics(assay(sce), 
                                      subsets = list(mt = is.mito))

if (args$filter_type == "manual") {
  qc <- metadata(sce)$qc_thresholds
  mt_percent <- rna.qc.metrics$subsets$mt * 100
  keep <- rna.qc.metrics$detected >= qc[qc$metric == "nFeature", "min"] &
    rna.qc.metrics$detected <= qc[qc$metric == "nFeature", "max"] &
    mt_percent < qc[qc$metric == "percent.mt", "max"] &
    rna.qc.metrics$sum <= qc[qc$metric == "nCount", "max"]
} else if (args$filter_type == "scrapper-auto") {
  require(DelayedArray)
  rna.qc.thresholds <- suggestRnaQcThresholds(rna.qc.metrics)
  keep <- filterRnaQcMetrics(rna.qc.thresholds, rna.qc.metrics)
}

output_file <- file.path(args$output_dir, paste0(args$name, "_cellids.txt.gz"))
writeLines(colnames(sce)[keep], output_file)

