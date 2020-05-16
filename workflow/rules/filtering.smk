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
