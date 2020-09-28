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

rule mqc_peaks_count:
    input:
        peaks="results/macs2_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak", # or narrowPeak if no --broad option
        header="../workflow/header/peaks_count_header.txt"
    output:
        "results/macs2_callpeak/peaks_count/{sample}-{control}.{peak}.peaks_count.tsv"
    log:
        "logs/macs2_callpeak/peaks_count/{sample}-{control}.{peak}.peaks_count.log"
    conda:
        "../envs/gawk.yaml"
    shell:
        "cat {input.peaks} | wc -l | gawk -v OFS='\t' '{{print \"{wildcards.sample}-{wildcards.control}_{wildcards.peak}_peaks\", $1}}' | cat {input.header} - > {output} 2> {log}"

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

rule create_mqc_frip_score:
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

rule create_igv_peaks:
    input:
        "results/macs2_callpeak/{sample}-{control}.{peak}_peaks.{peak}Peak"
    output:
        "results/IGV/macs2_callpeak/{peak}/merged_library.{sample}-{control}.{peak}_peaks.igv.txt"
    log:
        "logs/igv/create_igv_peaks/merged_library.{sample}-{control}.{peak}_peaks.log"
    shell:
        " find {input} -type f -name '*_peaks.{wildcards.peak}Peak' -exec echo -e 'results/IGV/macs2_callpeak/{wildcards.peak}/\"{{}}\"\t0,0,178' \; > {output} 2> {log}"
