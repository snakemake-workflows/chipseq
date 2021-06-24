`Base distribution by cycle plot
<https://gatk.broadinstitute.org/hc/en-us/articles/360042477312-CollectBaseDistributionByCycle-Picard->`_ **(Picard)** is
used as quality control for alignment-level and shows the nucleotide distribution per cycle of the bam files after
filtering, sorting, merging and removing orphans. For any cycle within reads the relative proportions of nucleotides
should reflect the AT:CG content. For all nucleotides flattish lines would be expected and any spikes would suggest a
systematic sequencing error. For more information about `collected Picard metrics
<https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectMultipleMetrics-Picard->`_ please
see `documentation <https://broadinstitute.github.io/picard/>`_.
