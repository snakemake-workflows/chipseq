`Plot of base distribution per sequencing cycle
<https://gatk.broadinstitute.org/hc/en-us/articles/360042477312-CollectBaseDistributionByCycle-Picard->`_.
This **Picard** tool shows the nucleotide distribution per sequencing cycle of the bam files after filtering, sorting, merging and removing orphans.
For any sequencing cycle within reads, the relative proportions of nucleotides should reflect the AT:CG content.
For all nucleotides, flat lines would be expected and any spikes suggest a systematic sequencing error.
For more information about `collected Picard metrics
<https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectMultipleMetrics-Picard->`_ please
see the `documentation <https://broadinstitute.github.io/picard/>`_.
