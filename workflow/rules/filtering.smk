rule samtools_view:
    input:
         get_samtools_view_input
    output:
        temp("results/sam-view/{sample}.bam")
    params:
        # TODO: move to config.yaml, including the explanatory comments below
        "-b -F 0x004 -G 0x009 -f 0x001 {} {}".format(get_blacklist_option(), get_blacklist_filter())
        # if duplicates should be removed in this filtering, add "-F 0x0400" to the params
        # if for each read, you only want to retain a single (best) mapping, add "-q 1" to params
        # if you would like to restrict analysis to certain regions (e.g. excluding other "blacklisted" regions),
        # please provide a respective bed file via "-L path/to/regions.bed"
    log:
        "logs/samtools-view/{sample}.log"
    wrapper:
        "0.64.0/bio/samtools/view"

rule bamtools_filter_json:
    input:
        "results/sam-view/{sample}.bam"
    output:
        temp("results/filtered/{sample}.bam")
    params:
          # filters mismatches in all reads and filters pe-reads within a size range given in json-file
        json="../config/bamtools_filtering_rules.json"
    log:
        "logs/filtered/{sample}.log"
    wrapper:
        "0.64.0/bio/bamtools/filter_json"

#TODO for later: customize and substitute rm_orphan_pe_bam.py with some existing tool
rule orphan_remove:
    input:
        "results/filtered/{sample}.bam"
    output:
        bam=temp("results/orphan_rm/{sample}.bam"),
        qc="results/orphan_rm/{sample}_bampe_rm_orphan.log"
    params:
        "--only_fr_pairs"
    log:
        "logs/orphan_remove/{sample}.log"
    conda:
        "../envs/pysam.yaml"
    shell:
        " ../workflow/scripts/rm_orphan_pe_bam.py {input} {output.bam} {params} 2> {log}"

rule samtools_sort:
    input:
        "results/orphan_rm/{sample}.bam"
    output:
        "results/orphan_rm_sorted/{sample}.bam"
    params:
        ""
    log:
        "logs/samtools-sort/{sample}.log"
    threads:  # Samtools takes additional threads through its option -@
        8
    wrapper:
        "0.64.0/bio/samtools/sort"
