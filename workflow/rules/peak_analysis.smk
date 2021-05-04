rule plot_fingerprint:
    input:
        bam_files=["results/filtered/{sample}.sorted.bam", "results/filtered/{control}.sorted.bam"],
        bam_idx=["results/filtered/{sample}.sorted.bam.bai", "results/filtered/{control}.sorted.bam.bai"],
        jsd_sample="results/filtered/{control}.sorted.bam",
        stats=expand("results/{step}/{{sample}}.sorted.{step}.stats.txt",
            step="bamtools_filtered" if config["single_end"]
            else "orph_rm_pe")
    output:  #ToDo: add description to report caption
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
        lambda w, input:
            "{se_option}{fragment_size}".format(
                se_option="--extendReads " if config["single_end"] else "",
                # Estimated fragment size used to extend single-end reads
                fragment_size=
                "$(grep ^SN {stats} | "
                "cut -f 2- | "
                "grep -m1 'average length:' | "
                "awk '{{print $NF}}') ".format(
                    stats=input.stats)
                if config["single_end"] else ""
            )
    threads:
        8
    wrapper:
        "0.66.0/bio/deeptools/plotfingerprint"

rule macs2_callpeak_broad:
    input:
        treatment="results/filtered/{sample}.sorted.bam",
        control="results/filtered/{control}.sorted.bam"
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
        "--broad-cutoff 0.1 -f {bam_format} {gsize} -B --SPMR --keep-dup all {pvalue} {qvalue}".format(
            gsize=get_gsize(),
            pvalue="-p {}".format(config["params"]["callpeak"]["p-value"]) if config["params"]["callpeak"][
                "p-value"] else "",
            qvalue="-p {}".format(config["params"]["callpeak"]["q-value"]) if config["params"]["callpeak"][
                "q-value"] else "",
            bam_format="BAM" if config["single_end"] else "BAMPE")
    wrapper:
        "0.68.0/bio/macs2/callpeak"

rule macs2_callpeak_narrow:
    input:
        treatment="results/filtered/{sample}.sorted.bam",
        control="results/filtered/{control}.sorted.bam"
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
        "-f {bam_format} {gsize} -B --SPMR --keep-dup all {pvalue} {qvalue}".format(
            gsize=get_gsize(),
            pvalue="-p {}".format(config["params"]["callpeak"]["p-value"]) if config["params"]["callpeak"][
                "p-value"] else "",
            qvalue="-p {}".format(config["params"]["callpeak"]["q-value"]) if config["params"]["callpeak"][
                "q-value"] else "",
            bam_format="BAM" if config["single_end"] else "BAMPE")
    wrapper:
        "0.68.0/bio/macs2/callpeak"

rule peaks_count:
    input:
        peaks="results/macs2_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak"
    output:
        "results/macs2_callpeak/peaks_count/{sample}-{control}.{peak}.peaks_count.tsv"
    log:
        "logs/macs2_callpeak/peaks_count/{sample}-{control}.{peak}.peaks_count.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "cat {input.peaks} | "
        " wc -l | "
        " gawk -v OFS='\t' '{{print \"{wildcards.sample}-{wildcards.control}_{wildcards.peak}_peaks\", $1}}' "
        " > {output} 2> {log}"

rule sm_report_peaks_count_plot:
    input:
        get_peaks_count_plot_input()
    output:
        report("results/macs2_callpeak/plots/plot_{peak}_peaks_count.pdf", caption="../report/plot_peaks_count_macs2.rst", category="CallPeaks")
    log:
        "logs/macs2_callpeak/plot_{peak}_peaks_count.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_peaks_count_macs2.R"

rule bedtools_intersect:
    input:
        left="results/filtered/{sample}.sorted.bam",
        right="results/macs2_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak"
    output:
        pipe("results/bedtools_intersect/{sample}-{control}.{peak}.intersected.bed")
    params:
        extra="-bed -c -f 0.20"
    log:
        "logs/bedtools/intersect/{sample}-{control}.{peak}.intersected.log"
    wrapper:
        "0.66.0/bio/bedtools/intersect"

rule frip_score:
    input:
        intersect="results/bedtools_intersect/{sample}-{control}.{peak}.intersected.bed",
        flagstats=expand("results/{step}/{{sample}}.sorted.{step}.flagstat", step= "bamtools_filtered" if config["single_end"]
        else "orph_rm_pe")
    output:
        "results/bedtools_intersect/{sample}-{control}.{peak}.peaks_frip.tsv"
    log:
        "logs/bedtools/intersect/{sample}-{control}.{peak}.peaks_frip.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "grep 'mapped (' {input.flagstats} | "
        " gawk -v a=$(gawk -F '\t' '{{sum += $NF}} END {{print sum}}' < {input.intersect}) "
        " -v OFS='\t' "
        " '{{print \"{wildcards.sample}-{wildcards.control}_{wildcards.peak}_peaks\", a/$1}}' "
        " > {output} 2> {log}"

rule sm_rep_frip_score:
    input:
        get_frip_score_input()
    output:
        report("results/macs2_callpeak/plots/plot_{peak}_peaks_frip_score.pdf", caption="../report/plot_frip_score_macs2_bedtools.rst", category="CallPeaks")
    log:
        "logs/bedtools/intersect/plot_{peak}_peaks_frip_score.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_frip_score.R"

rule create_igv_peaks:
    input:
        "results/macs2_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak"
    output:
        "results/IGV/macs2_callpeak-{peak}/merged_library.{sample}-{control}.{peak}_peaks.igv.txt"
    log:
        "logs/igv/create_igv_peaks/merged_library.{sample}-{control}.{peak}_peaks.log"
    shell:
        " find {input} -type f -name '*_peaks.{wildcards.peak}Peak' -exec echo -e 'results/IGV/macs2_callpeak/{wildcards.peak}/\"{{}}\"\t0,0,178' \; > {output} 2> {log}"

rule homer_annotatepeaks:
    input:
        peaks="results/macs2_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak",
        genome="resources/ref/genome.fasta",
        gtf="resources/ref/annotation.gtf"
    output:
        annotations="results/homer/annotate_peaks/{sample}-{control}.{peak}_peaks.annotatePeaks.txt"
    threads:
        2
    params:
        mode="",
        extra="-gid"
    log:
        "logs/homer/annotate_peaks/{sample}-{control}.{peak}.log"
    wrapper:
        "0.68.0/bio/homer/annotatePeaks"

rule plot_macs_qc:
    input:
        get_macs2_peaks()
    output:  #ToDo: add description to report caption
        summmary="results/macs2_callpeak/plots/plot_{peak}_peaks_macs2_summary.txt",
        plot=report("results/macs2_callpeak/plots/plot_{peak}_peaks_macs2.pdf", caption="../report/plot_macs2_qc.rst", category="CallPeaks")
    params:
        input = lambda wc, input: ','.join(input),
        sample_control_combinations = ','.join(get_sample_control_peak_combinations_list())
    log:
        "logs/macs2_callpeak/plot_{peak}_peaks_macs2.log"
    conda:
        "../envs/plot_macs_annot.yaml"
    shell:
        "Rscript ../workflow/scripts/plot_macs_qc.R -i {params.input} -s {params.sample_control_combinations}  -o {output.plot} -p {output.summmary} 2> {log}"

rule plot_homer_annotatepeaks:
    input:
        get_plot_homer_annotatepeaks_input()
    output:  #ToDo: add description to report caption
        summmary="results/homer/plots/plot_{peak}_annotatepeaks_summary.txt",
        plot=report("results/homer/plots/plot_{peak}_annotatepeaks.pdf", caption="../report/plot_annotatepeaks_homer.rst", category="CallPeaks")
    params:
        input = lambda wc, input: ','.join(input),
        sample_control_combinations = ','.join(get_sample_control_peak_combinations_list())
    log:
        "logs/homer/plot_{peak}_annotatepeaks.log"
    conda:
        "../envs/plot_macs_annot.yaml"
    shell:
        "Rscript ../workflow/scripts/plot_homer_annotatepeaks.R -i {params.input} -s {params.sample_control_combinations}  -o {output.plot} -p {output.summmary} 2> {log}"

rule plot_sum_annotatepeaks:
    input:
        "results/homer/plots/plot_{peak}_annotatepeaks_summary.txt"
    output:
        report("results/homer/plots/plot_{peak}_annotatepeaks_summary.pdf", caption="../report/plot_annotatepeaks_summary_homer.rst", category="CallPeaks")
    log:
        "logs/homer/plot_{peak}_annotatepeaks_summary.log"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_annotatepeaks_summary_homer.R"
