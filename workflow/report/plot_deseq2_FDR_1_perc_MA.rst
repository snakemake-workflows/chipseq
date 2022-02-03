The `MA plot <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#ma-plot>`_ **(FDR 0.01)**
displays the effect size (log2 fold changes) of differences in feature counts
between the group {{snakemake.wildcards["group_1"]}}
and the group {{snakemake.wildcards["group_2"]}}
(both treated with the {{snakemake.wildcards["antibody"]}} antibody)
as determined by DESeq2 (based on the mean of normalized counts).
The results of this plot have been filtered for a false discovery rate (FDR) of ``0.01``.
For more information about DESeq2 please see the
`documentation <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>`_.
