rule preseq_lc_extrap:
    input:
        "results/sam-view/{sample}.bam"
    output:
        "results/preseq/{sample}.lc_extrap"
    params:
        "-v"   #optional parameters
    log:
        "logs/preseq/{sample}.log"
    wrapper:
        "0.61.0/bio/preseq/lc_extrap"

rule collect_multiple_metrics:
    input:
         bam="results/orphan_rm_sorted/{sample}.bam",
         ref=config["resources"]["ref"]["genome"]
    output:
        # Through the output file extensions the different tools for the metrics can be selected
        # so that it is not necessary to specify them under params with the "PROGRAM" option.
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/collectmultiplemetrics.html.
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
        # This parameter (default 3 GB) can be used to limit the total resources a pipeline is allowed to use, see:
        #     https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources
        mem_gb=3
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

rule genomecov:
    input:
        "results/orphan_rm_sorted/{sample}.bam",
        "results/orphan_rm_sorted/{sample}.orphan_rm_sorted.flagstat"
    output:
        pipe("results/bed_graph/{sample}.bedgraph")
    log:
        "logs/bed_graph/{sample}.log"
    params: #-fs option was not used because there are no single end reads any more
        "-bg -pc -scale $(grep 'mapped (' results/orphan_rm_sorted/{sample}.orphan_rm_sorted.flagstat | awk '{print 1000000/$1}')"
    conda:
        "../envs/temp_genomecov.yaml"
    script:
        "../scripts/temp_bedtools_genomecoveragebed.py"
    #TODO: Add wrapper and remove script and env in workflow. Remove script and conda statements in this rule
    # wrapper:
    #     "xxxx/bio/bedtools/genomecov"

rule sort_genomecov:
    input:
        "results/bed_graph/{sample}.bedgraph"
    output:
        "results/bed_graph/{sample}.sorted.bedgraph"
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"

rule chromosome_size:
    input:
        config["resources"]["ref"]["index"]
    output:
        "ngs-test-data/ref/genome.chrom.sizes"
    shell:
        "cut -f 1,2 {input} > {output}"

rule bedGraphToBigWig:
    input:
        bedGraph="results/bed_graph/{sample}.sorted.bedgraph",
        chromsizes="ngs-test-data/ref/genome.chrom.sizes"
    output:
        "results/big_wig/{sample}.bigWig"
    params:
        ""
    wrapper:
        "0.61.0/bio/ucsc/bedGraphToBigWig"
