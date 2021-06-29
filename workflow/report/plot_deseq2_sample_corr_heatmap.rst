Correlation `heatmap plots <https://cran.r-project.org/web/packages/pheatmap/pheatmap.pdf>`_ shows heatmaps of the
`rlog <https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf#Rfn.rlog>`_ or
`vst <https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf#Rfn.vst>`_ transformed counts from
the DESeq2 analysis. The results of this plot represent the comparison of the
{{snakemake.wildcards["group_1"]}} versus {{snakemake.wildcards["group_2"]}} groups treated with the
{{snakemake.wildcards["antibody"]}} antibody. For more information about DESeq2 please
see `documentation <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>`_.
