rule cutadapt_pe:
    input:
        get_fastqs
    output:
        fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.pe.qc.txt",
        log="logs/cutadapt/{sample}-{unit}.pe.log"
    params:
        adapters = config["params"]["cutadapt-pe"],
        others = config["params"]["cutadapt-others"]
    log:
        "logs/cutadapt/{sample}-{unit}.pe.log"
    wrapper:
        "0.52.0/bio/cutadapt/pe"

rule cutadapt_se:
    input:
        get_fastqs
    output:
        fastq="results/trimmed/{sample}-{unit}.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.se.qc.txt",
        log="logs/cutadapt/{sample}-{unit}.se.log"
    params:
        "{} {}".format(
            config["params"]["cutadapt-se"],
            config["params"]["cutadapt-others"]
            )
    log:
        "logs/cutadapt/{sample}-{unit}.se.log"
    wrapper:
        "0.52.0/bio/cutadapt/se"
