#!/usr/bin/env Rscript
# Export DESeq2 VST to CSV for the morphology-to-RNA Python model.

args <- commandArgs(trailingOnly = TRUE)
vst_rds <- if (length(args) >= 1) args[[1]] else "/ems/elsc-labs/habib-n/segev.munitz/051223_BMP4_deseq2/data/vst_matrix.rds"
coldata_rds <- if (length(args) >= 2) args[[2]] else "/ems/elsc-labs/habib-n/segev.munitz/051223_BMP4_deseq2/data/colData.rds"
out_vst <- if (length(args) >= 3) args[[3]] else "data/vst_matrix.csv"
out_col <- if (length(args) >= 4) args[[4]] else "data/colData.csv"

vst <- readRDS(vst_rds)
if (inherits(vst, "DESeqTransform")) {
  mat <- as.matrix(SummarizedExperiment::assay(vst))
} else if (is.matrix(vst) || is.data.frame(vst)) {
  mat <- as.matrix(vst)
} else {
  stop("Unrecognized VST object class: ", paste(class(vst), collapse = ", "))
}

dir.create(dirname(out_vst), recursive = TRUE, showWarnings = FALSE)
write.csv(mat, out_vst, quote = TRUE, row.names = TRUE)

if (file.exists(coldata_rds)) {
  cd <- as.data.frame(readRDS(coldata_rds))
  cd$sample_id <- rownames(cd)
  write.csv(cd, out_col, quote = TRUE, row.names = FALSE)
}

message("Wrote ", out_vst, "  genes=", nrow(mat), " samples=", ncol(mat))
