rule get_genome:
    output:
        "resources/ref/genome.fasta"
    log:
        "logs/ref/get-genome.log"
    params:
        species=config["resources"]["ref"]["species"],
        datatype="dna",
        build=config["resources"]["ref"]["build"],
        release=config["resources"]["ref"]["release"],
        chromosome=config["resources"]["ref"]["chromosome"]
    cache: True
    wrapper:
        "0.67.0/bio/reference/ensembl-sequence"

rule get_annotation:
    output:
        "resources/ref/annotation.gtf"
    params:
        species=config["resources"]["ref"]["species"],
        release=config["resources"]["ref"]["release"],
        build=config["resources"]["ref"]["build"],
        fmt="gtf",
        flavor=""
    log:
        "logs/ref/get_annotation.log"
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "0.64.0/bio/reference/ensembl-annotation"

rule gtf2bed:
    input:
        "resources/ref/annotation.gtf"
    output:
        "resources/ref/genome.bed"
    log:
        "logs/ref/gtf2bed.log"
    conda:
        "../envs/perl.yaml"
    shell:
        "../workflow/scripts/gtf2bed {input} > {output} 2> {log}"

rule genome_faidx:
    input:
        "resources/ref/genome.fasta"
    output:
        "resources/ref/genome.fasta.fai"
    log:
        "logs/ref/genome-faidx.log"
    cache: True
    wrapper:
        "0.64.0/bio/samtools/faidx"

rule bwa_index:
    input:
        "resources/ref/genome.fasta"
    output:
        multiext("resources/ref/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa/bwa_index.log"
    cache: True
    params:
        algorithm="bwtsw"
    wrapper:
        "0.64.0/bio/bwa/index"

rule chromosome_size:
    input:
        "resources/ref/genome.fasta.fai"
    output:
        "resources/ref/genome.chrom.sizes"
    log:
        "logs/ref/chromosome_size.log"
    shell:
        "cut -f 1,2 {input} > {output} 2> {log}"

rule bedtools_sort_blacklist:
    input:
        in_file="../workflow/blacklists/{chrom}{blacklist}",
        sort_file="resources/ref/genome.chrom.sizes"
    output:
        "resources/ref/sorted_{chrom}{blacklist}"
    params:
        extra="-g"
    log:
        "logs/ref/sorted_{chrom}{blacklist}.log"
    # ToDo change to wrapper when released
    conda:
        "../envs/temp_bedtools_sort.yaml"
    script:
        "../scripts/temp_bedtools_sort.py"
    # wrapper:
    #     "xxx/bio/bedtools/sort"

rule bedtools_complement_blacklist:
    input:
        in_file="resources/ref/sorted_{chrom}{blacklist}",
        genome="resources/ref/genome.chrom.sizes"
    output:
        "resources/ref/sorted_complement_{chrom}{blacklist}"
    params:
        extra=""
    log:
        "logs/ref/sorted_complement_{chrom}{blacklist}.log"
        # ToDo change to wrapper when released
    conda:
        "../envs/temp_bedtools_sort.yaml"
    script:
        "../scripts/temp_bedtools_complement.py"
    # wrapper:
    #     "master/bio/bedtools/complement"

rule without_blacklists:
    input:
        "resources/ref/genome.chrom.sizes"
    output:
        "resources/ref/genome.chrom.sizes.filter"
    log:
        "logs/ref/genome.chrom.sizes.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "gawk '{{print \$1, \"0\" , \$2}}' OFS='\t' {input} > {output} 2> {log}"
