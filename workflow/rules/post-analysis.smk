rule preseq_lc_extrap:
    input:
        "results/sam-view/{sample}.bam"
    output:
        "results/preseq/{sample}.lc_extrap"
    params:
        "-bam -v"   #optional parameters
    log:
        "logs/preseq/{sample}.log"
    conda:
        "../envs/temp_preseq.yaml"
    script:
        "../scripts/temp_preseq_lc_extrap.py"
 #TODO: add wrapper and remove script, env and cript and conda statement in this rule
    # wrapper:
    #     "xxxx/bio/preseq/lc_extrap"

rule alignment_summary:
    input:
        ref=config["resources"]["ref"]["genome"],
        bam="results/orphan_rm_sorted/{sample}.bam"
    output:
        "results//alignment-summary/{sample}.summary.txt"
    log:
        "logs/picard/alignment-summary/{sample}.log"
    params:
        # optional parameters
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE"
    wrapper:
        "0.60.1/bio/picard/collectalignmentsummarymetrics"

