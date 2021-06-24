import os

from snakemake import rules

rule bedtools_merge_broad:
    input:
        get_macs2_peaks()
    output:
        "results/bedtools/merged/{antibody}.consensus_broad-peaks.txt"
    params:
        extra="-c {} -o {}".format( ','.join(map(str, list( range(2,10) ) ) ),
                                       ','.join( ["collapse"] * 8))
    log:
        "logs/bedtools/merged/{antibody}.consensus_peaks.log"
    wrapper:
        "0.66.0/bio/bedtools/merge"

rule bedtools_merge_narrow:
    input:
        get_macs2_peaks()
    output:
        "results/bedtools/merged/{antibody}.consensus_narrow-peaks.txt"
    params:
        extra="-c {} -o {}".format( ','.join(map(str, list( range(2,11) ) ) ),
                                       ','.join( ["collapse"] * 9))
    log:
        "logs/bedtools/merged/{antibody}.consensus_peaks.log"
    wrapper:
        "0.66.0/bio/bedtools/merge"

rule macs2_merged_expand:
    input:
        "results/bedtools/merged/{antibody}.consensus_{peak}-peaks.txt"
    output:
        bool_txt="results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.txt",
        bool_intersect="results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.intersect.txt"
    params:
        sample_control_peak=get_sample_control_peak_combinations_list(),
        narrow_param="--is_narrow_peak" if config["params"]["peak-analysis"] == "narrow" else "",
        min_reps_consensus=config["params"]["min-reps-consensus"]
    log:
        "logs/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.log"
    script:
        "../scripts/macs2_merged_expand.py"

rule create_consensus_bed:
    input:
        "results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.txt"
    output:
        "results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.bed"
    conda:
        "../envs/gawk.yaml"
    log:
        "logs/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.bed.log"
    shell:
        "gawk -v FS='\t' -v OFS='\t' 'FNR  > 1 {{ print $1, $2, $3, $4 \"0\", \"+\"}}' {input} > {output} 2> {log}"

rule create_consensus_saf:
    input:
        "results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.txt"
    output:
        "results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.saf"
    conda:
        "../envs/gawk.yaml"
    log:
        "logs/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.bed.log"
    shell:
        "$(echo -e 'GeneID\tChr\tStart\tEnd\tStrand' > {output} && "
        " gawk -v FS='\t' -v OFS='\t' 'FNR > 1 {{ print $4, $1, $2, $3,  \" + \" }}' {input} >> {output}) "
        " 2> {log}"

rule plot_peak_intersect:
    input:
        "results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.intersect.txt"
    output:
       report(
           "results/macs2_merged_expand/plots/{antibody}.consensus_{peak}-peaks.boolean.intersect.plot.pdf",
           caption="../report/plot_consensus_peak_intersect.rst",
           category="ConsensusPeak")
    conda:
        "../envs/consensus_plot.yaml"
    log:
        "logs/macs2_merged_expand/plots/{antibody}.consensus_{peak}-peaks.boolean.intersect.plot.log"
    shell:
        "Rscript ../workflow/scripts/plot_peak_intersect.R -i {input} -o {output} 2> {log}"

rule create_consensus_igv:
    input:
        "results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.bed"
    output:
        "results/IGV/consensus/merged_library.{antibody}.consensus_{peak}-peaks.igv.txt"
    log:
        "logs/igv/consensus/merged_library.{antibody}.consensus_{peak}-peaks.igv.log"
    shell:
        "find {input} -type f -name '*.consensus_{wildcards.peak}-peaks.boolean.bed' -exec "
        "echo -e '{{}}\t0,0,0' \; > {output} 2> {log}"

rule homer_consensus_annotatepeaks:
    input:
        peaks="results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.bed",
        genome="resources/ref/genome.fasta",
        gtf="resources/ref/annotation.gtf"
    output:
        annotations="results/homer/annotate_consensus_peaks/{antibody}.consensus_{peak}-peaks.annotatePeaks.txt"
    threads:
        2
    resources:
        mem_mb=4000
    params:
        mode="",
        extra="-gid"
    log:
        "logs/homer/annotate_consensus_peaks/{antibody}.consensus_{peak}-peaks.annotatePeaks.log"
    wrapper:
        "0.68.0/bio/homer/annotatePeaks"

rule trim_homer_consensus_annotatepeaks:
    input:
        "results/homer/annotate_consensus_peaks/{antibody}.consensus_{peak}-peaks.annotatePeaks.txt"
    output:
        temp("results/homer/annotate_consensus_peaks/{antibody}.consensus_{peak}-peaks.annotatePeaks.trimmed.txt")
    log:
        "logs/homer/annotate_consensus_peaks/trimmed/{antibody}.consensus_{peak}-peaks.annotatePeaks.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "cut -f2- {input} | gawk 'NR==1; NR > 1 {{print $0 | \"sort -T '.' -k1,1 -k2,2n\"}}' | cut -f6- > {output}"

rule merge_bool_and_annotatepeaks:
    input:
        trim="results/homer/annotate_consensus_peaks/{antibody}.consensus_{peak}-peaks.annotatePeaks.trimmed.txt",
        bool="results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.txt"
    output:
        "results/homer/annotate_consensus_peaks/{antibody}.consensus_{peak}-peaks.boolean.annotatePeaks.txt"
    log:
        "logs/homer/annotate_consensus_peaks/{antibody}.consensus_{peak}-peaks.boolean.annotatePeaks.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "paste {input.bool} {input.trim} > {output}"

rule feature_counts:
    input:
        sam=lambda wc: expand("results/filtered/{sample}.sorted.bam",
            sample=get_samples_of_antibody(wc.antibody)),
        annotation="results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.saf"
    output:
        multiext("results/feature_counts/{antibody}.consensus_{peak}-peaks",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts")
    threads:
        2
    params:
        extra="-F SAF -O --fracOverlap 0.2{pe_param}".format(pe_param="" if config["single_end"] else " -p --donotsort")
    log:
        "logs/feature_counts/{antibody}.consensus_{peak}-peaks.featureCounts.log"
    wrapper:
        "0.73.0/bio/subread/featurecounts"

rule featurecounts_modified_colnames:
    input:
        featurecounts="results/feature_counts/{antibody}.consensus_{peak}-peaks.featureCounts",
        bam=expand("results/filtered/{sample}.sorted.bam", sample=samples.index),
        samples_file=config["samples"]
    output:
        "results/feature_counts/{antibody}.consensus_{peak}-peaks_modified.featureCounts"
    params:
        ""
    log:
        "logs/feature_counts/{antibody}.consensus_{peak}-peaks_modified.featureCounts.log"
    script:
        "../scripts/col_mod_featurecounts.py"

rule featurecounts_deseq2:
    input:
        "results/feature_counts/{antibody}.consensus_{peak}-peaks_modified.featureCounts"
    output:
        dds="results/deseq2/dss_rld/{antibody}.consensus_{peak}-peaks.dds.rld.RData",
        plot_pca=report("results/deseq2/plots/{antibody}.consensus_{peak}-peaks.pca_plot.pdf",
            caption = "../report/plot_deseq2_pca.rst", category = "DESeq2"),
        plot_heatmap=report("results/deseq2/plots/{antibody}.consensus_{peak}-peaks.heatmap_plot.pdf",
            caption = "../report/plot_deseq2_heatmap.rst", category = "DESeq2"),
        pca_data="results/deseq2/pca_vals/{antibody}.consensus_{peak}-peaks.pca.vals.txt",
        dist_data="results/deseq2/dists/{antibody}.consensus_{peak}-peaks.sample.dists.txt",
        size_factors_rdata="results/deseq2/sizeFactors/{antibody}.consensus_{peak}-peaks.sizeFactors.RData",
        size_factors_res="results/deseq2/sizeFactors/{antibody}.consensus_{peak}-peaks.sizeFactors.sizeFactor.txt",
        results="results/deseq2/results/{antibody}.consensus_{peak}-peaks.deseq2_results.txt",
        # pairwise comparisons of samples across the groups from a particular antibody
        FDR_1_perc_res=directory("results/deseq2/FDR/results/FDR_0.01_{antibody}.consensus_{peak}-peaks"),
        FDR_5_perc_res=directory("results/deseq2/FDR/results/FDR_0.05_{antibody}.consensus_{peak}-peaks"),
        FDR_1_perc_bed=directory("results/deseq2/FDR/bed_files/FDR_0.01_{antibody}.consensus_{peak}-peaks"),
        FDR_5_perc_bed=directory("results/deseq2/FDR/bed_files/FDR_0.05_{antibody}.consensus_{peak}-peaks"),
        igv_FDR_5_bed="results/IGV/consensus/merged_library.{antibody}.consensus_{peak}-peaks.deseq2.FDR_0.05.igv.txt",
        plot_FDR_1_perc_MA=report(
            directory("results/deseq2/comparison_plots/MA_plots/FDR_0.01_{antibody}consensus_{peak}-peaks"),
            patterns=["{antibody}.{{group_1_vs_group_2}}.MA-plot_FDR_0.01.pdf"],
            caption = "../report/plot_deseq2_FDR_1_perc_MA.rst",
            category = "DESeq2"),
        plot_FDR_5_perc_MA=report(
            directory("results/deseq2/comparison_plots/MA_plots/FDR_0.05_{antibody}consensus_{peak}-peaks"),
            patterns=["{antibody}.{{group_1_vs_group_2}.MA-plot_FDR_0.05.pdf"],
            caption = "../report/plot_deseq2_FDR_5_perc_MA.rst",
            category = "DESeq2"),
        plot_FDR_1_perc_volcano=report(
            directory("results/deseq2/comparison_plots/volcano_plots/FDR_0.01_{antibody}consensus_{peak}-peaks"),
            patterns=["{antibody}.{{group_1_vs_group_2}}.volcano-plot_FDR_0.01.pdf"],
            caption = "../report/plot_deseq2_FDR_1_perc_volcano.rst",
            category = "DESeq2"),
        plot_FDR_5_perc_volcano=report(
            directory("results/deseq2/comparison_plots/volcano_plots/FDR_0.05_{antibody}consensus_{peak}-peaks"),
            patterns=["{antibody}.{{group_1_vs_group_2}}.volcano-plot_FDR_0.05.pdf"],
            caption = "../report/plot_deseq2_FDR_5_perc_volcano.rst",
            category = "DESeq2"),
        plot_sample_corr_heatmap=report(
            directory("results/deseq2/comparison_plots/correlation_heatmaps_{antibody}consensus_{peak}-peaks"),
            patterns=["{antibody}.{{group_1_vs_group_2}}.correlation_heatmap.pdf"],
            caption = "../report/plot_deseq2_sample_corr_heatmap.rst",
            category = "DESeq2"),
        plot_scatter=report(
            directory("results/deseq2/comparison_plots/scatter_plots_{antibody}consensus_{peak}-peaks"),
            patterns=["{antibody}.{{group_1_vs_group_2}}.scatter_plots.pdf"],
            caption = "../report/plot_deseq2_scatter.rst",
            category = "DESeq2")
    threads:
        2
    params:
        vst = config["params"]["deseq2"]["vst"],
        antibody = lambda w: w.antibody
    log:
        "logs/deseq2/{antibody}.consensus_{peak}-peaks.featureCounts.log"
    conda:
        "../envs/featurecounts_deseq2.yaml"
    script:
        "../scripts/featurecounts_deseq2.R"
