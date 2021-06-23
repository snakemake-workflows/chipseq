rule fastqc:
    input:
        get_individual_fastq
    output:
        html="results/qc/fastqc/{sample}.{unit}.{read}.html",
        zip="results/qc/fastqc/{sample}.{unit}.{read}_fastqc.zip"
    params:
        ""
    log:
        "logs/fastqc/{sample}.{unit}.{read}.log"
    threads: 6
    wrapper:
        "0.72.0/bio/fastqc"

rule multiqc:
    input:
        get_multiqc_input
    output:
        report("results/qc/multiqc/multiqc.html", category="MultiQC report")
    log:
        "logs/multiqc.log"
    wrapper:
        "0.64.0/bio/multiqc"
