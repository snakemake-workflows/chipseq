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
