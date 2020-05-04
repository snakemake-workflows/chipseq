rule cutadapt_pe:
    input:
        get_fastqs
    output:
        fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.pe.qc.txt"
    params:
        adapters = "-a AGAGCACACGTCTGAACTCCAGTCAC -g AGATCGGAAGAGCACACGT -A AGAGCACACGTCTGAACTCCAGTCAC -G AGATCGGAAGAGCACACGT",
        others = "--minimum-length 1 -q 20"
        # "{}".format(config["params"]["cutadapt-pe"])
    log:
        "results/logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.52.0/bio/cutadapt/pe"

rule cutadapt_se:
    input:
        get_fastqs
    output:
        fastq="results/trimmed/{sample}-{unit}.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.se.qc.txt"
    params:
        "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 20"
        # "{}".format(config["params"]["cutadapt-se"])
    log:
        "results/logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.52.0/bio/cutadapt/se"
