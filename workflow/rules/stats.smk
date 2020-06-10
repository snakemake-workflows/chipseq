rule samtools_flagstat:
    input:
        "results/{path_and_file_name}.bam"
    output:
        "results/{path_and_file_name}.flagstat"
    wrapper:
        "0.60.0/bio/samtools/flagstat"

rule samtools_idxstats:
    input:
        bam = "results/{path_and_file_name}.bam",
        idx = "results/{path_and_file_name}.bam.bai"
    output:
        "results/{path_and_file_name}.idxstats"
    log:
        "logs/{path_and_file_name}.log"
    wrapper:
        "0.60.0/bio/samtools/idxstats"

rule samtools_stats:
    input:
        "results/{path_and_file_name}.bam"
    output:
        "results/{path_and_file_name}.stats.txt"
    params:
        ""
    log:
        "logs/{path_and_file_name}.log"
    wrapper:
        "0.60.0/bio/samtools/stats"
