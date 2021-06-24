rule create_igv_input_file:
    input:
        get_igv_input()
    output:
        temp("results/IGV/igv_files.txt")
    log:
        "logs/igv/merged_igv_files.log"
    shell:
        "cat {input} > {output} 2> {log}"

# igv session that can be started directly from the generated files of workflow
rule igv_files_to_session:
    input:
        igv="results/IGV/igv_files.txt",
        fasta="resources/ref/genome.fasta"
    output:
        "results/IGV/igv_session.xml"
    params:
        "--path_prefix '../../'"
    log:
        "logs/igv/igv_session_to_file.log"
    shell:
        " ../workflow/scripts/igv_files_to_session.py {output} {input.igv} ../../{input.fasta} {params} 2> {log}"

# remove paths for igv session to download from report.zip
rule igv_report_cleanup:
    input:
        "results/IGV/igv_files.txt"
    output:
        temp("results/IGV/report_igv_files.txt")
    log:
        "logs/igv/igv_files.log"
    script:
        "../scripts/igv_report_cleanup.py"

rule igv_files_to_report:
    input:
        igv="results/IGV/report_igv_files.txt",
        fasta="resources/ref/genome.fasta"
    output:
        temp("results/IGV/report_igv_session.xml")
    params:
        ""
    log:
        "logs/igv/igv_files_to_report.log"
    shell:
        " ../workflow/scripts/igv_files_to_session.py {output} {input.igv} $(basename {input.fasta}) {params} 2> {log}"

rule collect_igv_report_session_files:
    input:
        igv_data=get_files_for_igv(),
        deseq2_files=directory(expand("results/deseq2/FDR/bed_files/FDR_0.05_{antibody}.consensus_{peak}-peaks",
        antibody=get_unique_antibodies(), peak=config["params"]["peak-analysis"])),
        fasta="resources/ref/genome.fasta",
        igv_session="results/IGV/report_igv_session.xml"
    output:
        temp(directory("results/IGV/report_igv_session"))
    log:
        "logs/igv/collect_igv_report_session_files.log"
    shell:
        "mkdir -p {output}; cp {input.igv_data} {output}/; "
        "cp -R {input.deseq2_files}/* {output}/; cp {input.fasta} {output}/; "
        "cp {input.igv_session} {output}/"

# igv session that can be downloaded from generated report
rule zip_igv_report_session:
    input:
        rules.collect_igv_report_session_files.output
    output:
        report("results/IGV/report_igv_session.zip", caption = "../report/igv_session.rst", category="IGV session")
    log:
        "logs/igv/collect_igv_report_session_files.log"
    conda:
        "../envs/gzip.yaml"
    shell:
        "cd $(dirname {input}); zip $(basename {output}) $(basename {input})/*"
