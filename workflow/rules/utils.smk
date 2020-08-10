rule samtools_index:
    input:
        "results/{step}/{samples_units}.bam"
    output:
        "results/{step}/{samples_units}.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.64.0/bio/samtools/index"
