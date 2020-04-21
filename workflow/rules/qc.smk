rule fastqc:
    input:
        expand("{path}/{{read}}.fq", path=config["reads_dir"])
    output:
        html="results/qc/fastqc/reports/{read}.fq.html",
        zip="results/qc/fastqc/{read}.fq_fastqc.zip"
    params: ""
    log:
        "logs/fastqc/{read}.log"
    wrapper:
        "0.51.2/bio/fastqc"
