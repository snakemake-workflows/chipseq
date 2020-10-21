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
