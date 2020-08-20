rule fastqc:
    input:
        get_individual_fastq
    output:
        html="results/qc/fastqc/{sample}.{unit}.{read}.html",
        zip="results/qc/fastqc/{sample}.{unit}.{read}_fastqc.zip"
    log:
        "logs/fastqc/{sample}.{unit}.{read}.log"
    wrapper:
        "0.64.0/bio/fastqc"

# ToDo: add wrapper again and remove temporary script and env after wrapper release with matplotlib dependency
rule multiqc:
    input:
        get_multiqc_input
    output:
         "results/qc/multiqc/multiqc.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "0.64.0/bio/multiqc"
