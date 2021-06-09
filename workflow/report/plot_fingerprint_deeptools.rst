**deepTools `plotFingerprint <https://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html>`_** is a
quality control tool, which determines how well the signal in the ChIP-seq sample can be differentiated from the
background distribution of reads in the control sample. plotFingerprint randomly samples genome regions of a specified
length and sums the per-base coverage that overlap with those regions. These values are then sorted according to their
rank and the cumulative sum of read counts is plotted. With a perfect uniform distribution of reads along the genome and
infinite sequencing coverage a straight diagonal line is shown in the plot. For more information on the interpretation of
the plots please see `here <https://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html#what-the-plots-tell-you>`_.
For more information about plotFingerprint metrics see
`here <https://deeptools.readthedocs.io/en/latest/content/feature/plotFingerprint_QC_metrics.html>`_.
