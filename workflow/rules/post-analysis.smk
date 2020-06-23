rule preseq_lc_extrap:
    input:
        "results/sam-view/{sample}.bam"
    output:
        "results/preseq/{sample}.lc_extrap"
    params:
        "-v"   #optional parameters
    log:
        "logs/preseq/{sample}.log"
    conda:
        "../envs/temp_preseq.yaml"
    script:
        "../scripts/temp_preseq_lc_extrap.py"
 #TODO: Add wrapper and remove script and env in workflow. Remove script and conda statements in this rule
    # wrapper:
    #     "xxxx/bio/preseq/lc_extrap"

rule collect_multiple_metrics:
    input:
         bam="results/orphan_rm_sorted/{sample}.bam",
         ref=config["resources"]["ref"]["genome"]
    output:
        multiext("{path}{sample}",
                 ".alignment_summary_metrics",
                 ".base_distribution_by_cycle_metrics",
                 ".base_distribution_by_cycle.pdf",
                 ".insert_size_metrics",
                 ".insert_size_histogram.pdf",
                 ".quality_by_cycle_metrics",
                 ".quality_by_cycle.pdf",
                 ".quality_distribution_metrics",
                 ".quality_distribution.pdf"
                 )
    resources:
        mem_gb=3
    wildcard_constraints:
        # the common path for all metrics must be defined here
        path="results/qc/multiple_metrics/"
    log:
        "logs/picard/{path}{sample}.log"
    params:
        # optional parameters
        "VALIDATION_STRINGENCY=LENIENT "
    conda:
        "../envs/temp_collectmultiplemetrics.yaml"
    script:
        "../scripts/temp_picard_collectmultiplemetrics.py"
    #TODO: Add wrapper and remove script and env in workflow. Remove script and conda statements in this rule
    # wrapper:
    #     "xxxx/bio/picard/collectmultiplemetrics"
