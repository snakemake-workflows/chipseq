rule bwa_index:
    input:
        config["resources"]["ref"]["genome"]
    output:
        multiext(config["resources"]["ref"]["genome"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa/bwa_index.log"
    params:
        algorithm="bwtsw"
    wrapper:
        "0.55.1/bio/bwa/index"

rule bwa_mem:
    input:
        reads = get_map_reads_input,
        idx = rules.bwa_index.output
    output:
        temp("results/mapped/{sample}-{unit}.bam")
    log:
        "logs/bwa/bwa_mem/{sample}-{unit}.log"
    params:
        index= lambda w, input: os.path.splitext(input.idx[0])[0],
        extra= get_read_group,
        sort="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "0.52.0/bio/bwa/mem"

rule mark_duplicates:
    input:
        rules.bwa_mem.output
    output:
        bam=temp("results/mapped/dedup/{sample}-{unit}.bam"),
        metrics="results/mapped/dedup/{sample}-{unit}.metrics.txt"
    log:
        "logs/picard/dedup_{sample}-{unit}.log"
    params:
        "REMOVE_DUPLICATES=true"
    wrapper:
        "0.55.1/bio/picard/markduplicates"

rule merge_bams:
    input:
        get_dedup_bams
    output:
        temp("results/mapped/merged_with_dups.bam")
    log:
        "logs/picard/mergebamfiles.log"
    params:
        "VALIDATION_STRINGENCY=LENIENT"
    wrapper:
        "0.55.1/bio/picard/mergesamfiles"

rule mark_merged_duplicates:
    input:
        "results/mapped/merged_with_dups.bam"
    output:
        bam="results/mapped/merged.bam",
        metrics="results/mapped/dedup/merged.metrics.txt"
    log:
        "logs/picard/dedup_merged.log"
    params:
        "REMOVE_DUPLICATES=true"
    wrapper:
        "0.55.1/bio/picard/markduplicates"
