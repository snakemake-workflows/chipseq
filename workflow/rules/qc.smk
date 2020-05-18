rule fastqc:
    input:
        get_individual_fastq
    output:
        html="results/qc/fastqc/{sample}.{unit}.{read}.html",
        zip="results/qc/fastqc/{sample}.{unit}.{read}_fastqc.zip"
    log:
        "logs/fastqc/{sample}.{unit}.{read}.log"
    wrapper:
        "0.51.2/bio/fastqc"

rule multiqc:
    input:
        get_multiqc_input
    output:
         "results/qc/multiqc/multiqc.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "0.51.3/bio/multiqc"
