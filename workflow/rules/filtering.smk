rule samtools_view:
    input:
        "results/merged/dedup/{sample}.bam"
    output:
        temp("results/filtered/sam-view/{sample}.bam")  ###TODO pipe to bamtools_filter_json
    params:
        "-b -F 0x004 -G 0x009 -f 0x001"
        # if duplicates should be removed in this filtering, add "-F 0x0400" to the params
        # if for each read, you only want to retain a single (best) mapping, add "-q 1" to params
        # if you would like to restrict analysis to certain regions (e.g. excluding other "blacklisted" regions),
        # please provide a respective bed file via "-L path/to/regions.bed"
    wrapper:
        "0.60.0/bio/samtools/view"

rule bamtools_filter_json:
    input:
        "results/filtered/sam-view/{sample}.bam"
    output:
        temp("results/filtered/bamtools/{sample}.bam")
    params:
          # filters mismatches in all reads and filters pe-reads within a size range given in json-file
        json="../config/bamtools_filtering_rules.json"
    log:
        "logs/filtered/bamtools/{sample}.log"
    wrapper:
        "0.60.0/bio/bamtools/filter_json"

rule samtools_index_filtered:
    input:
        "results/filtered/bamtools/{sample}.bam"
    output:
        "results/filtered/bamtools/{sample}.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.60.0/bio/samtools/index"

rule samtools_flagstat_filtered:
    input:
        "results/filtered/bamtools/{sample}.bam"
    output:
        "results/filtered/flagstat/filtered_{sample}.flagstat"
    wrapper:
        "0.60.0/bio/samtools/flagstat"

rule idxstats_filtered:
    input:
        bam = "results/filtered/bamtools/{sample}.bam",
        idx = "results/filtered/bamtools/{sample}.bam.bai"
    output:
        "results/filtered/idxstats/filtered_{sample}.idxstats"
    log:
        "logs/samtools/idxstats/filtered_{sample}.log"
    wrapper:
        "0.60.0/bio/samtools/idxstats"

rule samtools_stats_filtered:
    input:
        "results/filtered/bamtools/{sample}.bam"
    output:
        "results/filtered/stats/filtered_{sample}.stats.txt"
    params:
        ""
    log:
        "logs/samtools/stats/filtered_{sample}.log"
    wrapper:
        "0.60.0/bio/samtools/stats"

rule orphan_remove:
    input:
        "results/filtered/bamtools/{sample}.bam"
    output:
        bam=temp("results/filtered/orphan_rm/{sample}.bam"),
        qc="results/filtered/orphan_rm/{sample}_bampe_rm_orphan.log"
    conda:
        "../envs/pysam.yaml"
    shell:
        " ../workflow/scripts/rm_orphan_pe_bam.py {input} {output.bam} --only_fr_pairs "

rule samtools_sort:
    input:
        "results/filtered/orphan_rm/{sample}.bam"
    output:
        "results/filtered/sorted/{sample}.sorted.bam"
    params:
        ""
    threads:  # Samtools takes additional threads through its option -@
        workflow.cores * 0.75
    wrapper:
        "0.60.0/bio/samtools/sort"

rule samtools_index_rm_orphan:
    input:
        "results/filtered/sorted/{sample}.sorted.bam"
    output:
        "results/filtered/sorted/{sample}.sorted.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.60.0/bio/samtools/index"

rule samtools_flagstat_rm_orphan:
    input:
        "results/filtered/sorted/{sample}.sorted.bam"
    output:
        "results/filtered/flagstat/sorted_{sample}.flagstat"
    wrapper:
        "0.60.0/bio/samtools/flagstat"

rule idxstats_rm_orphan:
    input:
        bam = "results/filtered/sorted/{sample}.sorted.bam",
        idx = "results/filtered/sorted/{sample}.sorted.bam.bai"
    output:
        "results/filtered/idxstats/sorted_{sample}.idxstats"
    log:
        "logs/samtools/idxstats/sorted_{sample}.log"
    wrapper:
        "0.60.0/bio/samtools/idxstats"

rule samtools_stats_rm_orphan:
    input:
        "results/filtered/sorted/{sample}.sorted.bam"
    output:
        "results/filtered/stats/sorted_{sample}.stats.txt"
    params:
        ""
    log:
        "logs/samtools/stats/sorted_{sample}.log"
    wrapper:
        "0.60.0/bio/samtools/stats"
