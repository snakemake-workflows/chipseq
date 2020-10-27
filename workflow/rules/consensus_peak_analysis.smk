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
