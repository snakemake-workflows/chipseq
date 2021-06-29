**Volcano plot (FDR 0.05)** shows the significance (adjusted p-value) versus the log2 fold changes of the
DESeq2 analysis results.
The results of this plot are filtered on a false discovery rate (FDR) threshold of 0.05 and represent the comparison of the
{{snakemake.wildcards["group_1"]}} versus {{snakemake.wildcards["group_2"]}} groups for the
{{snakemake.wildcards["antibody"]}} antibody. For more information about DESeq2 please see
`documentation <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>`_.
