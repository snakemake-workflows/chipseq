SAMPLE_CONTROL=get_sample_control_combinations()

rule plot_fingerprint:
    input:
        bam_files=["results/orphan_rm_sorted/{sample}.bam", "results/orphan_rm_sorted/{control}.bam"],
        bam_idx=["results/orphan_rm_sorted/{sample}.bam.bai", "results/orphan_rm_sorted/{control}.bam.bai"],
        jsd_sample="results/orphan_rm_sorted/{control}.bam"
    output:
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotfingerprint.html.
        fingerprint="results/deeptools/{sample}-{control}.plot_fingerprint.pdf",
        counts="results/deeptools/{sample}-{control}.fingerprint_counts.txt",
        qc_metrics="results/deeptools/{sample}-{control}.fingerprint_qcmetrics.txt"
    log:
        "logs/deeptools/plot_fingerprint.{sample}-{control}.log"
    params:
        "--labels {sample} {control}",
        "--skipZeros ",
        "--numberOfSamples 500000 ", # ToDo: to config?
    threads:
        8
    wrapper:
        "0.66.0/bio/deeptools/plotfingerprint"

rule macs2_callpeak_broad:
    input:
        treatment="results/orphan_rm_sorted/{sample}.bam",
        control="results/orphan_rm_sorted/{control}.bam"
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("results/macs2_callpeak/{sample}-{control}.callpeak",
                 "_peaks.xls",
                 # these output extensions internally set the --bdg or -B option:
                 "_treat_pileup.bdg",
                 "_control_lambda.bdg",
                 # these output extensions internally set the --broad option:
                 "_peaks.broadPeak",
                 "_peaks.gappedPeak"
                 )
    log:
        "logs/macs2/callpeak.{sample}-{control}.log"
    params: # ToDo: make --broad option (and other options?) selectable in config
        "--broad --broad-cutoff 0.1 -f BAMPE -g hs --SPMR --qvalue 0.05 --keep-dup all"
        # ToDo: Update wrapper to check for " --broad$" or " --broad " instead of only "--broad" (line 47),
        #  then "--broad" in params can be removed here in params
    wrapper:
        "0.66.0/bio/macs2/callpeak"

# rule mqc_peaks_count:
#     input:
#         peaks="results/macs2_callpeak/{sample}-{control}.callpeak_peaks.broadPeak", # or narrowPeak if no --broad option
#         header="../workflow/header/peaks_count_header.txt"
#     output:
#         "results/macs2_callpeak/stats/{sample}-{control}.peaks_count.tsv"
#     log:
#         "logs/macs2_callpeak/{sample}-{control}.peaks_count.log"
#     conda:
#         "../envs/gawk.yaml"
#     shell:
#         "( cat {input.peaks} | wc -l | gawk -v OFS='\t' '{{print \"{wildcards.sample}-{wildcards.control}\", $1}}' | cat {input.header} - > {output} 2>&1 ) >{log}"
#

rule mqc_peaks_count:
    input:
        peaks=expand("results/macs2_callpeak/{sam_contr}.callpeak_peaks.broadPeak", sam_contr=SAMPLE_CONTROL), # or narrowPeak if no --broad option
        header="../workflow/header/peaks_count_header.txt"
    output:
        "results/macs2_callpeak/stats/peaks_count.tsv"
    log:
        "logs/macs2_callpeak/peaks_count.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "( cat {input.header} > {output}; for i in {input.peaks}; do cat $i | wc -l | gawk -v OFS='\t' '{{print $i, $1}}' >> {output}; done 2>&1 ) >{log}"


rule bedtools_intersect:
    input:
        left="results/orphan_rm_sorted/{sample}.bam",
        right="results/macs2_callpeak/{sample}-{control}.callpeak_peaks.broadPeak"
    output:
        "results/intersect/{sample}-{control}.intersected.bed"
    params:
        extra="-bed -c -f 0.20"
    log:
        "logs/intersect/{sample}-{control}.intersected.log"
    wrapper:
        "0.66.0/bio/bedtools/intersect"
