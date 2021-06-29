The `MA plot <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#ma-plot>`_ **(FDR 0.01)**
displays the log2 fold changes versus the mean of normalized counts of the
DESeq2 analysis results.
The results of this plot are filtered on a false discovery rate (FDR) threshold of 0.01 and represent the comparison of the
{{snakemake.wildcards["group_1"]}} versus {{snakemake.wildcards["group_2"]}} groups treated with the
{{snakemake.wildcards["antibody"]}} antibody. For more information about DESeq2 please see
`documentation <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>`_.
