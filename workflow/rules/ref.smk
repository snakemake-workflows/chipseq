# TODO: reinstate the wrapper once it can handle downloading only a specific chromosome
# rule get_genome:
#     output:
#         "resources/ref/genome.fasta"
#     log:
#         "logs/ref/get-genome.log"
#     params:
#         species=config["resources"]["ref"]["species"],
#         datatype="dna",
#         build=config["resources"]["ref"]["build"],
#         release=config["resources"]["ref"]["release"]
#     cache: True
#     wrapper:
#         "0.64.0/bio/reference/ensembl-sequence"

rule get_genome_test_data:
    output:
        "resources/ref/genome.fasta"
    log:
        "logs/ref/get-genome.log"
    params:
        species=config["resources"]["ref"]["species"],
        spec_up=config["resources"]["ref"]["species"].capitalize(),
        datatype="dna",
        build=config["resources"]["ref"]["build"],
        suffix="chromosome.21",
        release=config["resources"]["ref"]["release"]
    cache: True
    conda:
        "../envs/curl.yaml"
    shell:
        "curl -L ftp://ftp.ensembl.org/pub/release-{params.release}/fasta/{params.species}/{params.datatype}/{params.spec_up}.{params.build}.{params.datatype}.{params.suffix}.fa.gz |gzip -d > {output} 2> {log}"

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
