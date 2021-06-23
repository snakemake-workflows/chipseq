rule create_igv_input_file:
    input:
        get_igv_input()
    output:
        # "results/IGV/merged_igv_files.txt"
        "results/IGV/igv_files.txt"
    log:
        "logs/igv/merged_igv_files.log"
    shell:
        "cat {input} > {output} 2> {log}"

# remove paths for igv session download from report.zip
# rule igv_report_cleanup:
#     input:
#         "results/IGV/merged_igv_files.txt"
#     output:
#         "results/IGV/igv_files.txt"
#     log:
#         "logs/igv/igv_files.log"
#     script:
#         " ../scripts/igv_report_cleanup.py"

rule igv_files_to_session:
    input:
        # igv_data=get_igv_data,
        igv="results/IGV/igv_files.txt",
        fasta="resources/ref/genome.fasta"
    output:
        report("results/IGV/igv_session.xml", category = "igv session")
    params:
        "--path_prefix '../../'"
    log:
        "logs/igv/igv_session_to_file.log"
    shell:
        " ../workflow/scripts/igv_files_to_session.py {output} {input.igv} ../../{input.fasta} {params} 2> {log}"
