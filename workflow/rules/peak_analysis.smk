rule plot_fingerprint:
    input:
        bam_files=get_sample_control_input,
        bam_idx=get_sample_control_idx
    output:
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotfingerprint.html.
        fingerprint="results/deeptools/plot_fingerprint.{sample}-{control}.pdf",
        counts="results/deeptools/fingerprint_counts.{sample}-{control}.txt",
        qc_metrics="results/deeptools/fingerprint_qcmetrics.{sample}-{control}.txt"
    log:
        "logs/deeptools/plot_fingerprint.{sample}-{control}.log"
    params:
        "--labels {sample} {control}",
        "--outQualityMetrics results/deeptools/fingerprint_qcmetrics.{sample}-{control}.txt",
        "--skipZeros ",
        "--numberOfProcessors 8",
        "--numberOfSamples 100 ", # default: 500000
        "--JSDsample results/orphan_rm_sorted/{control}.bam"
    threads:
        8
    wrapper:
        "0.64.0/bio/deeptools/plotfingerprint"
