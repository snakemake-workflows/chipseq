# bam_files=expand("results/orphan_rm_sorted/{sample}.bam", sample=samples.index) # for JSDsample in ToDo 3
rule plot_fingerprint:
    input:
         bam_files=expand("results/orphan_rm_sorted/{sample}.bam", sample=samples.index),
         bam_idx=expand("results/orphan_rm_sorted/{sample}.bam.bai", sample=samples.index)
    output:
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotfingerprint.html.
        fingerprint="results/deeptools/plot_fingerprint.pdf",
        counts="results/deeptools/fingerprint_counts.txt",
        qc_metrics="results/deeptools/fingerprint_qcmetrics.txt"
    log:
        "logs/deeptools/plot_fingerprint.log"
    params:
        "--labels "+" ".join(samples.index),
        "--outQualityMetrics results/deeptools/fingerprint_qcmetrics.txt ", # ToDo 1: change Wrapper & meta: additional output for qc-metrics
        "--skipZeros ",
        "--numberOfProcessors 8", # ToDo 2: change Wrapper: add threads, e.g. threads = "" if snakemake.threads <= 1 else " -@ {} ".format(snakemake.threads - 1)
        "--numberOfSamples 100 ", # default: 500000
    # " --JSDsample "+" --JSDsample ".join(bam_files), # ToDo 3: integrate into the wrapper?
    threads:
        8  # integrate in --numberOfProcessors -> ToDo 2
    wrapper:
        "0.64.0/bio/deeptools/plotfingerprint"

rule testing:
    input:
        get_sample_and_control
    output:
        "results/test.txt"
    shell:
        "echo {input} > {output}"
