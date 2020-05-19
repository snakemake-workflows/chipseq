rule samtools_view:
    input:
        "results/merged/dedup/{sample}.bam"
    output:
        "results/filtered/sam-view/{sample}.bam"
    params:
        "-b -F 0x004 -G 0x009 -f 0x001"
        # if duplicates should be removed in this filtering, add "-F 0x0400" to the params
        # if for each read, you only want to retain a single (best) mapping, add "-q 1" to params
        # if you would like to restrict analysis to certain regions (e.g. excluding other "blacklisted" regions),
        # please provide a respective bed file via "-L path/to/regions.bed"
    wrapper:
        "0.57.0/bio/samtools/view"



###TODO wrapper for bamtools sorting -> output: "results/filtered/sorted/{sample}.sorted.bam"

# rule samtools_index:
#     input:
#         "results/filtered/sorted/{sample}.sorted.bam"
#     output:
#         "results/filtered/sorted/{sample}.sorted.bam.bai"
#     params:
#         "" # optional params string
#     wrapper:
#         "0.57.0/bio/samtools/index"
#
# rule samtools_flagstat:
#     input:
#         "results/filtered/sorted/{sample}.sorted.bam"
#     output:
#         "results/filtered/sorted/flagstat/{sample}.bam.flagstat"
#     wrapper:
#         "0.57.0/bio/samtools/flagstat"

# rule idxstats:
#     input:
#         bam = "results/filtered/sorted/{sample}.sorted.bam",
#         idx = "results/filtered/sorted/{sample}.sorted.bam.bai"
#     output:
#         "results/filtered/sorted/idxstats/{sample}.bam.idxstats"
#     log:
#         "logs/samtools/idxstats/{sample}.log"
#     wrapper:
#         "0.57.0/bio/samtools/idxstats"

# rule samtools_stats:
#     input:
#         "results/filtered/sorted/{sample}.sorted.bam"
#     output:
#         "results/filtered/sorted/stats/{sample}.txt"
#     params:
#         extra="",                       # Optional: extra arguments.
#         region="xx:1000000-2000000"      # Optional: region string.
#     log:
#         "logs/samtools/stats/{sample}.log"
#     wrapper:
#         "0.57.0/bio/samtools/stats"
