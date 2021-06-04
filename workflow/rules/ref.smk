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

# SRA-download
rule sra_get_fastq_pe:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "resources/ref/sra-pe-reads/{accession}.1.fastq",
        "resources/ref/sra-pe-reads/{accession}.2.fastq"
    params:
        extra=""
    threads: 6
    log:
        "logs/ref/sra-pe-reads/{accession}.log"
    wrapper:
        "0.72.0/bio/sra-tools/fasterq-dump"

rule sra_get_fastq_se:
    output:
        "resources/ref/sra-se-reads/{accession}.fastq"
    params:
        extra=""
    threads: 6
    log:
        "logs/ref/sra-pe-reads/{accession}.log"
    wrapper:
        "0.72.0/bio/sra-tools/fasterq-dump"

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
        genome="resources/ref/genome.fasta.fai"
    output:
        "resources/ref/genome.chrom.sizes"
    log:
        "logs/ref/chromosome_size.log"
    shell:
        "cut -f 1,2 {input.genome} > {output} 2> {log}"

rule generate_igenomes:
    output:
        "resources/ref/igenomes.yaml"
    params:
        igenomes_release = config["resources"]["ref"]["igenomes_release"]
    log:
        "logs/ref/igenomes.log"
    conda:
        ""
    script:
        "../scripts/generate_igenomes.py"

rule generate_igenomes_blacklist:
    input:
        "resources/ref/igenomes.yaml"
    output:
        blacklist_path="resources/ref/blacklist.bed"
    params:
        build = config["resources"]["ref"]["build"],
        chromosome = config["resources"]["ref"]["chromosome"],
        blacklist = config["resources"]["ref"]["blacklist"]
    log:
        "logs/ref/blacklist.log"
    conda:
        ""
    script:
        "../scripts/generate_blacklist.py"

rule bedtools_sort_blacklist:
    input:
        in_file="resources/ref/blacklist.bed"
    output:
        "resources/ref/blacklist.sorted"
    params:
        extra=""
    log:
        "logs/ref/blacklist.sorted.log"
    wrapper:
        "0.68.0/bio/bedtools/sort"

rule bedtools_complement_blacklist:
    input:
        in_file="resources/ref/blacklist.sorted",
        genome="resources/ref/genome.chrom.sizes"
    output:
        "resources/ref/blacklist.sorted.complement"
    params:
        extra=""
    log:
        "logs/ref/blacklist.sorted.complement.log"
    wrapper:
        "0.68.0/bio/bedtools/complement"

checkpoint get_gsize:
    input:
        "resources/ref/igenomes.yaml"
    output:
        "resources/ref/gsize.txt"
    params:
        extra=config["resources"]["ref"]["macs-gsize"],
        build=config["resources"]["ref"]["build"]
    log:
        "logs/ref/gsize.log"
    conda:
        ""
    script:
        "../scripts/get_gsize.py"
