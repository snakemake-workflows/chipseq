`PCA plot <https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf#Rfn.plotPCA>`_
**(Principal Component Analysis)** show the variance of the
{% if snakemake.params.vst == "true" %}vst{% else %}rlog{% endif %} transformed feature counts
from DESeq2 with regard to the used {{snakemake.wildcards["antibody"]}} antibody.
There is more detailed `information on the count transformation in the docs <https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations>`_.
For more information about DESeq2 in general, please also see the
`documentation <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>`_.
