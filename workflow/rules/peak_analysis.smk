rule plot_fingerprint:
    input:
        bam_files=["results/orphan_rm_sorted/{sample}.bam", "results/orphan_rm_sorted/{control}.bam"],
        bam_idx=["results/orphan_rm_sorted/{sample}.bam.bai", "results/orphan_rm_sorted/{control}.bam.bai"],
        jsd_sample="results/orphan_rm_sorted/{control}.bam"
    output:
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotfingerprint.html.
        fingerprint="results/deeptools/{sample}-{control}.plot_fingerprint.pdf",
        counts="results/deeptools/{sample}-{control}.fingerprint_counts.txt",
        qc_metrics="results/deeptools/{sample}-{control}.fingerprint_qcmetrics.txt"
    log:
        "logs/deeptools/plot_fingerprint.{sample}-{control}.log"
    params:
        "--labels {sample} {control}",
        "--skipZeros ",
        "--numberOfSamples 500000 ", # ToDo: to config?
    threads:
        8
    wrapper:
        "0.66.0/bio/deeptools/plotfingerprint"

rule macs2_callpeak_broad:
    input:
        treatment="results/orphan_rm_sorted/{sample}.bam",
        control="results/orphan_rm_sorted/{control}.bam"
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("results/macs2_callpeak/{sample}-{control}.callpeak",
                 "_peaks.xls",
                 # these output extensions internally set the --bdg or -B option:
                 "_treat_pileup.bdg",
                 "_control_lambda.bdg",
                 # these output extensions internally set the --broad option:
                 "_peaks.broadPeak",
                 "_peaks.gappedPeak"
                 )
    log:
        "logs/macs2/callpeak.{sample}-{control}.log"
    params: # ToDo: make --broad option (and other options?) selectable in config
        "--broad --broad-cutoff 0.1 -f BAMPE -g hs --SPMR --qvalue 0.05 --keep-dup all"
        # ToDo: Update wrapper to check for " --broad$" or " --broad " instead of only "--broad" (line 47),
        #  then "--broad" in params can be removed here in params
    wrapper:
        "0.66.0/bio/macs2/callpeak"
