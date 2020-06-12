rule samtools_flagstat:
    input:
        "results/{step}/{samples_units}.bam"
    output:
        "results/{step,[^./]+}/{samples_units}.{step}.flagstat"
    wrapper:
        "0.60.0/bio/samtools/flagstat"

rule samtools_idxstats:
    input:
        bam = "results/{step}/{samples_units}.bam",
        idx = "results/{step}/{samples_units}.bam.bai"
    output:
        "results/{step,[^./]+}/{samples_units}.{step}.idxstats"
    log:
        "logs/{step,[^./]+}/{samples_units}.{step}.log"
    wrapper:
        "0.60.0/bio/samtools/idxstats"

rule samtools_stats:
    input:
        "results/{step}/{samples_units}.bam"
    output:
        "results/{step,[^./]+}/{samples_units}.{step}.stats.txt"
    params:
        ""
    log:
        "logs/{step,[^./]+}/{samples_units}.{step}.log"
    wrapper:
        "0.60.0/bio/samtools/stats"
