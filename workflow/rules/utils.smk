rule samtools_index:
    input:
        "results/{path_and_file_name}.bam"
    output:
        "results/{path_and_file_name}.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.60.0/bio/samtools/index"
