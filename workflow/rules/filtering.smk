rule samtools_view:
    input:
         get_samtools_view_input
    output:
        temp("results/sam-view/{sample}.bam")
    params:
        # if duplicates should be removed in this filtering, add "-F 0x0400" to the params
        # if for each read, you only want to retain a single (best) mapping, add "-q 1" to params
        # if you would like to restrict analysis to certain regions (e.g. excluding other "blacklisted" regions),
        # the -L option is automatically activated if a path to a blacklist of the given genome exists in
        # ".test/config/igenomes.yaml" or has been entered there
        lambda wc, input: "-b -F 0x004 {pe_params} {bl}".format(
            pe_params="" if config["single_end"] else "-G 0x009 -f 0x001",
            bl="" if len(input) == 1 else "-L {}".format(list(input)[1])
        )
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
        json="../config/{}_bamtools_filtering_rules.json".format("se" if config["single_end"] else "pe")
    log:
        "logs/filtered/{sample}.log"
    wrapper:
        "0.64.0/bio/bamtools/filter_json"

rule samtools_sort_pe:
    input:
        "results/filtered/{sample}.bam"
    output:
        "results/filtered/{sample}.sorted.bam"
    params:
        ""
    log:
        "logs/samtools-sort/{sample}.sorted.log"
    threads:
        8
    wrapper:
        "0.64.0/bio/samtools/sort"

#TODO for later: customize and substitute rm_orphan_pe_bam.py with some existing tool
rule orphan_remove:
    input:
        "results/filtered/{sample}.sorted.bam"
        # expand("results/filtered/{{sample}}{infix}.bam",
        #    infix="" if config["single_end"] else ".sorted"
       # )
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
         expand("{path}/{{sample}}.bam",
                path="results/filtered" if config["single_end"] else "results/orphan_rm")
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
