#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
system(paste0("cp ", snakemake@input[["header"]], " ", snakemake@output[[1]]))
load(snakemake@input[["data"]])
write.table(crosscorr['cross.correlation'], file=snakemake@output[[1]], sep=',', quote=FALSE,
row.names=FALSE, col.names=FALSE, append=TRUE)
