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
        temp("results/merged/{sample}-{unit}.bam")
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

rule merge_bams:
    input:
        lambda w: expand("results/merged/{sample}-{unit}.bam",
            sample = w.sample,
            unit = units.loc[units['sample'] == w.sample].unit.to_list()
        )
    output:
        temp("results/merged/{sample}.bam")
    log:
        "logs/picard/mergebamfiles/{sample}.log"
    params:
        "VALIDATION_STRINGENCY=LENIENT"
    wrapper:
        "0.55.1/bio/picard/mergesamfiles"

rule mark_merged_duplicates:
    input:
        "results/merged/{sample}.bam"
    output:
        bam="results/merged/dedup/{sample}.bam",
        metrics="results/merged/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup_merged/{sample}.log"
    params:
        "REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT"
    wrapper:
        "0.55.1/bio/picard/markduplicates"
