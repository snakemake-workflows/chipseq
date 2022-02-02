The **Volcano plot (FDR 0.01)** shows the significance (adjusted p-value) versus the effect size (log2 fold changes) 
differences in feature counts between the group {{snakemake.wildcards["group_1"]}}
and the group {{snakemake.wildcards["group_2"]}} group
(both treated with the {{snakemake.wildcards["antibody"]}} antibody)
as determined by DESeq2 (based on the mean of normalized counts).
The results of this plot have been filtered for a false discovery rate (FDR) of ``0.05``.
For more information about DESeq2 please see the
`documentation <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>`_.
