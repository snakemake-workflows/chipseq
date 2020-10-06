rule plot_fingerprint:
    input:
        bam_files=["results/orphan_rm_sorted/{sample}.bam", "results/orphan_rm_sorted/{control}.bam"],
        bam_idx=["results/orphan_rm_sorted/{sample}.bam.bai", "results/orphan_rm_sorted/{control}.bam.bai"],
        jsd_sample="results/orphan_rm_sorted/{control}.bam"
    output:
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotfingerprint.html.
        fingerprint=report("results/deeptools/{sample}-{control}.plot_fingerprint.pdf", caption="../report/plot_fingerprint_deeptools.rst", category="QC"),
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
        multiext("results/macs2_callpeak/{sample}-{control}.broad",
                 "_peaks.xls",
                 # these output extensions internally set the --bdg or -B option:
                 "_treat_pileup.bdg",
                 "_control_lambda.bdg",
                 # these output extensions internally set the --broad option:
                 "_peaks.broadPeak",
                 "_peaks.gappedPeak"
                 )
    log:
        "logs/macs2/callpeak.{sample}-{control}.broad.log"
    params: # ToDo: move to config?
        "--broad --broad-cutoff 0.1 -f BAMPE -g hs --SPMR --qvalue 0.05 --keep-dup all"
        # ToDo: Update wrapper to check for " --broad$" or " --broad " instead of only "--broad" (line 47),
        #  then "--broad" in params can be removed here in params
    wrapper:
        "0.66.0/bio/macs2/callpeak"

rule macs2_callpeak_narrow:
    input:
        treatment="results/orphan_rm_sorted/{sample}.bam",
        control="results/orphan_rm_sorted/{control}.bam"
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("results/macs2_callpeak/{sample}-{control}.narrow",
                 "_peaks.xls",
                 # these output extensions internally set the --bdg or -B option:
                 "_treat_pileup.bdg",
                 "_control_lambda.bdg",
                 # these output extensions internally set the --broad option:
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
    log:
        "logs/macs2/callpeak.{sample}-{control}.narrow.log"
    params: # ToDo: move to config?
        "-f BAMPE -g hs --SPMR --qvalue 0.05 --keep-dup all"
    wrapper:
        "0.66.0/bio/macs2/callpeak"

rule mqc_peaks_count:  # not displayed in multiqc.html -> added to snakemake-report plot_peaks_count_macs2.rst, see next rule
    input:
        peaks="results/macs2_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak",
        header="../workflow/header/peaks_count_header.txt"
    output:
        "results/macs2_callpeak/peaks_count/{sample}-{control}.{peak}.peaks_count.tsv"
    log:
        "logs/macs2_callpeak/peaks_count/{sample}-{control}.{peak}.peaks_count.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "cat {input.peaks} | wc -l | gawk -v OFS='\t' '{{print \"{wildcards.sample}-{wildcards.control}_{wildcards.peak}_peaks\", $1}}' | cat {input.header} - > {output} 2> {log}"

rule sm_report_peaks_count:
    input:
        expand("results/macs2_callpeak/peaks_count/{sam_contr}.{peak}.peaks_count.tsv",
               sam_contr=get_sample_control_combinations(), peak=list(config["params"]["peak_analysis"].split()))
    output:
        report("results/macs2_callpeak/plots/plot_peaks_count.pdf", caption="../report/plot_peaks_count_macs2.rst", category="CallPeaks")
    log:
        "logs/macs2_callpeak/plot_peaks_count.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_peaks_count_macs2.R"

rule bedtools_intersect:
    input:
        left="results/orphan_rm_sorted/{sample}.bam",
        right="results/macs2_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak"
    output:
        pipe("results/intersect/{sample}-{control}.{peak}.intersected.bed")
    params:
        extra="-bed -c -f 0.20"
    log:
        "logs/intersect/{sample}-{control}.{peak}.intersected.log"
    wrapper:
        "0.66.0/bio/bedtools/intersect"

rule create_mqc_frip_score:  # not displayed in multiqc.html -> added to snakemake-report plot_frip_score_macs2_bedtools.rst, see next rule
    input:
        intersect="results/intersect/{sample}-{control}.{peak}.intersected.bed",
        flagstats="results/orphan_rm_sorted/{sample}.orphan_rm_sorted.flagstat",
        header="../workflow/header/frip_score_header.txt"
    output:
        "results/intersect/{sample}-{control}.{peak}.peaks_frip.tsv"
    log:
        "logs/intersect/{sample}-{control}.{peak}.peaks_frip.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "grep 'mapped (' {input.flagstats} | gawk -v a=$(gawk -F '\t' '{{sum += $NF}} END {{print sum}}' < {input.intersect}) -v OFS='\t' "
        "'{{print \"{wildcards.sample}-{wildcards.control}_{wildcards.peak}_peaks\", a/$1}}' | cat {input.header} - > {output} 2> {log}"

# rule sm_rep_frip_score:
#     input:
#         expand("results/intersect/{sam_contr}.{peak}.peaks_frip.tsv",
#                sam_contr=get_sample_control_combinations(), peak=config["params"]["peak_analysis"])
#     output:
#         report("results/intersect/plot_peaks_frip_score.pdf", caption="../report/plot_frip_score_macs2_bedtools.rst", category="CallPeaks")
#     log:
#         "logs/intersect/plot_peaks_frip_score.log"
#     conda:
#         "../envs/r_plots.yaml"
#     script:
#         "../scripts/plot_frip_score_macs2_bedtools.R"

rule create_igv_peaks:
    input:
        "results/macs2_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak"
    output:
        "results/IGV/macs2_callpeak/{peak}/merged_library.{sample}-{control}.{peak}_peaks.igv.txt"
    log:
        "logs/igv/create_igv_peaks/merged_library.{sample}-{control}.{peak}_peaks.log"
    shell:
        " find {input} -type f -name '*_peaks.{wildcards.peak}Peak' -exec echo -e 'results/IGV/macs2_callpeak/{wildcards.peak}/\"{{}}\"\t0,0,178' \; > {output} 2> {log}"

rule homer_annotatepeaks:
    input:
        peaks="results/macs2_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak",
        genome="resources/ref/genome.fasta",
        gtf="resources/ref/annotation.gtf",
    output:
        annotations="results/homer/annotate_peaks/{sample}-{control}.{peak}_peaks.annotatePeaks.txt"
    threads:
        2
    params:
        mode="",
        extra="-gid"
    log:
        "logs/homer/annotate_peaks/{sample}-{control}.{peak}.log"
    conda:
        "../envs/temp_annotatepeaks.yaml"
    script:
        "../scripts/temp_annotatepeaks_wrapper.py"
    #  #TODO: add wrapper and remove script, env and conda statement in this rule
    # wrapper:
    #     "xxx/bio/homer/annotatePeaks"

rule plot_macs_qc:
    input:
        expand("results/macs2_callpeak/{sam_contr}.{peak}_peaks.{peak}Peak",
               sam_contr=get_sample_control_combinations(), peak =config["params"]["peak_analysis"])
    output:
        summmary="results/macs2_callpeak/plots/macs2_plot_summary.txt",
        plot=report("results/macs2_callpeak/plots/macs2_plot.pdf", caption="../report/plot_macs2_qc.rst", category="CallPeaks")
    params:
        get_sample_control_peak_combinations()
    log:
        "logs/macs2_callpeak/macs2_plot.log"
    conda:
        "../envs/plot_macs_annot.yaml"
    shell:
        "Rscript ../workflow/scripts/plot_macs_qc.R -i $(echo {input} | sed 's/ /,/g') -s $(echo {params}| sed 's/ /,/g')  -o {output.plot} -p {output.summmary} 2> {log}"

rule plot_homer_annotatepeaks:
    input:
        expand("results/homer/annotate_peaks/{sam_contr}.{peak}_peaks.annotatePeaks.txt",
               sam_contr=get_sample_control_combinations(), peak =config["params"]["peak_analysis"])
    output:
        summmary="results/homer/plots/annotatepeaks_plot_summary.txt",
        plot=report("results/homer/plots/annotatepeaks_plot.pdf", caption="../report/plot_annotatepeaks_homer.rst", category="CallPeaks")
    params:
        get_sample_control_peak_combinations()
    log:
        "logs/homer/annotatepeaks_plot.log"
    conda:
        "../envs/plot_macs_annot.yaml"
    shell:
        "Rscript ../workflow/scripts/plot_homer_annotatepeaks.R -i $(echo {input} | sed 's/ /,/g') -s $(echo {params}| sed 's/ /,/g')  -o {output.plot} -p {output.summmary} 2> {log}"

rule create_mqc_plot_annotatepeaks: # not displayed in multiqc.html -> added to snakemake-report plot_annotatepeaks_summary_homer.rst, see next rule
    input:
        summary="results/homer/plots/annotatepeaks_plot_summary.txt",
        header="../workflow/header/peak_annotation_header.txt"
    output:
        "results/homer/plots/mqc_annotatepeaks_plot_summary.tsv"
    log:
        "logs/homer/mqc_annotatepeaks.log"
    shell:
        " cat {input.header} {input.summary} > {output} 2> {log}"

# rule sm_report_plot_sum_annotatepeaks:
#     input:
#         "results/homer/plots/mqc_annotatepeaks_plot_summary.tsv"
#     output:
#         report("results/homer/plots/plot_annotatepeaks_summary.pdf", caption="../report/plot_annotatepeaks_summary_homer.rst", category="CallPeaks")
#     log:
#         "logs/homer/plot_annotatepeaks_summary.log"
#     conda:
#         "../envs/r_plots.yaml"
#     script:
#         "../scripts/plot_annotatepeaks_summary_homer.R"
