rule plot_fingerprint:
    input:
        bam_files=["results/orphan_rm_sorted/{sample}.bam", "results/orphan_rm_sorted/{control}.bam"],
        bam_idx=["results/orphan_rm_sorted/{sample}.bam.bai", "results/orphan_rm_sorted/{control}.bam.bai"]
    output:
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotfingerprint.html.
        fingerprint="results/deeptools/plot_fingerprint.{sample}-{control}.pdf",
        counts="results/deeptools/fingerprint_counts.{sample}-{control}.txt",
        qc_metrics="results/deeptools/fingerprint_qcmetrics.{sample}-{control}.txt"
    log:
        "logs/deeptools/plot_fingerprint.{sample}-{control}.log"
    params:
        "--labels {sample} {control}",
        "--outQualityMetrics results/deeptools/fingerprint_qcmetrics.{sample}-{control}.txt", # ToDo: remove on wrapper update
        "--skipZeros ",
        "--numberOfProcessors 8", # ToDo: remove on wrapper update
        "--numberOfSamples 100 ", # default: 500000
        "--JSDsample results/orphan_rm_sorted/{control}.bam" # ToDo: move to input on wrapper update
    threads:
        8
    wrapper:
        "0.64.0/bio/deeptools/plotfingerprint"
