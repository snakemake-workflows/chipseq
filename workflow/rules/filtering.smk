rule samtools_view:
    input:
         get_samtools_view_input
    output:
        temp("results/sam-view/{sample}.bam")
    params:
        "{sv_params} {bl} {bl_filter}".format(sv_params=get_samtools_view_params(), bl=get_blacklist_option(), bl_filter=get_blacklist_filter())
    log:
        "logs/samtools-view/{sample}.log"
    wrapper:
        "0.64.0/bio/samtools/view"

rule bamtools_filter_json:
    input:
        "results/sam-view/{sample}.bam"
    output:
        temp("results/filtered/{{sample}}{}.bam".format(get_pe_prefix()))
    params:
          # filters mismatches in all reads and filters pe-reads within a size range given in json-file
        json="../config/{}_bamtools_filtering_rules.json".format(get_se_pe_prefix())
    log:
        "logs/filtered/{sample}.log"
    wrapper:
        "0.64.0/bio/bamtools/filter_json"

rule samtools_sort_pe:
    input:
        "results/filtered/{sample}-pe.bam"
    output:
        "results/filtered/{sample}-pe.bam"
    params:
        ""
    log:
        "logs/samtools-sort/{sample}-pe.log"
    threads:
        8
    wrapper:
        "0.64.0/bio/samtools/sort"

#TODO for later: customize and substitute rm_orphan_pe_bam.py with some existing tool
rule orphan_remove:
    input:
        "results/filtered/{{sample}}{}.bam".format(get_pe_prefix())
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
