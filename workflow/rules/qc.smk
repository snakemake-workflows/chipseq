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

rule bwa_mem:
    input:
        reads = get_fastqs
        # reads=["reads/{sample}.1.fastq", "reads/{sample}.2.fastq"]
    output:
        "results/mapped/{sample}-{unit}.bam"
    log:
        "logs/bwa_mem/{sample}-{unit}.log"
    params:
        index= "genome",  #config["resources"]["ref"]["genome"],
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        # sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "0.52.0/bio/bwa/mem"
