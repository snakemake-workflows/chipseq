rule bedtools_merge_broad:
    input:
        get_plot_macs_qc_input()
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
        get_plot_macs_qc_input()
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
        narrow_param=get_narrow_flag(),
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
       report("results/macs2_merged_expand/plots/{antibody}.consensus_{peak}-peaks.boolean.intersect.plot.pdf", caption="../report/plot_consensus_peak_intersect.rst", category="ConsensusPeak")
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
        "find {input} -type f -name '*.consensus_{wildcards.peak}-peaks.boolean.bed' -exec echo -e 'results/IGV/consensus/{wildcards.antibody}/\"{{}}\"\t0,0,0' \; > {output} 2> {log}"

rule homer_consensus_annotatepeaks:
    input:
        peaks="results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.bed",
        genome="resources/ref/genome.fasta",
        gtf="resources/ref/annotation.gtf"
    output:
        annotations="results/homer/annotate_consensus_peaks/{antibody}.consensus_{peak}-peaks.annotatePeaks.txt"
    threads:
        2
    params:
        mode="",
        extra="-gid"
    log:
        "logs/homer/annotate_consensus_peaks/{antibody}.consensus_{peak}-peaks.annotatePeaks.log"
    # ToDo change to wrapper when released
    conda:
        "../envs/temp_annotatepeaks.yaml"
    script:
        "../scripts/temp_annotatepeaks_wrapper.py"
    # wrapper:
    #     "xxx/bio/homer/annotatePeaks"

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
        sam=get_samples_of_antibody,
        annotation="results/macs2_merged_expand/{antibody}.consensus_{peak}-peaks.boolean.saf"
        # optional input:
        # fasta="",     # implicitly sets the -G flag
        # chr_names=""  # implicitly sets the -A flag
    output:
        multiext("results/feature_counts/{antibody}.consensus_{peak}-peaks",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts")
    threads:
        2
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-F SAF -O --fracOverlap 0.2"
    log:
        "logs/feature_counts/{antibody}.consensus_{peak}-peaks.featureCounts.log"
    conda:
        "../envs/temp_featurecounts.yaml"
    script:
        "../scripts/temp_featurecounts.py"
    # ToDo change to wrapper when released
    # wrapper:
    #     "0.67.0/bio/homer/annotatePeaks"
