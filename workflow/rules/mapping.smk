rule bwa_mem:
    input:
        reads = get_map_reads_input,
        idx = rules.bwa_index.output
    output:
        "results/mapped/{sample}-{unit}.bam"
    log:
        "logs/bwa/bwa_mem/{sample}-{unit}.log"
    params:
        index= lambda w, input: os.path.splitext(input.idx[0])[0],
        extra= get_read_group,
        sort="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "0.64.0/bio/bwa/mem"

rule merge_bams:
    input:
        lambda w: expand("results/mapped/{sample}-{unit}.bam",
            sample = w.sample,
            unit = units.loc[units['sample'] == w.sample].unit.to_list()
        )
    output:
        "results/merged/{sample}.bam"
    log:
        "logs/picard/mergebamfiles/{sample}.log"
    params:
        "VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate TMP_DIR=tmp"
    wrapper:
        "0.64.0/bio/picard/mergesamfiles"

rule mark_merged_duplicates:
    input:
        "results/merged/{sample}.bam"
    output:
        bam="results/picard_dedup/{sample}.bam",
        metrics="results/picard_dedup/{sample}.metrics.txt"
    log:
        "logs/picard/picard_dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp"
    wrapper:
        "0.64.0/bio/picard/markduplicates"
