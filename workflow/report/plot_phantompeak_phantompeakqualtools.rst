 **`Phantompeakqualtools plot <https://code.google.com/archive/p/phantompeakqualtools/>`_** shows strand-shift versus
cross-correlation and computes informative enrichment and quality measures for ChIP-seq data. It also calculates the
relative (RSC) and the normalized strand cross-correlation coefficient (NSC). Datasets with NSC values much less than
1.1 tend to have low signal to noise or few peaks. This may indicate poor quality or only few binding sites. RSC values
significantly lower than 1 (< 0.8) tend to have low signal to noise. The low scores can be due to several factors and
are often due to failed and poor quality ChIP or low read sequence quality. In section "ChIP-seq processing pipeline"
of MultiQC report are additional plots across all samples for RSC, NSC and strand-shift cross-correlation.
 For more information about Phantompeakqualtools please
see `documentation <https://github.com/kundajelab/phantompeakqualtools/blob/master/README.md>`_. For more information
about interpretation see published `article <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/>`_.
