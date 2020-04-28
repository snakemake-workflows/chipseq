rule samples_fq:
    input:
        get_individual_fastq
    output:
        "results/samples_fq/{sample}.{unit}.{read}.fq"
    shell:
        "cp {input} {output}"


rule fastqc:
    input:
        "results/samples_fq/{sample}.{unit}.{read}.fq"
    output:
        html="results/qc/fastqc/reports/{sample}.{unit}.{read}.fq.html",
        zip="results/qc/fastqc/zip-files/{sample}.{unit}.{read}.fq_fastqc.zip"
    params: ""
    log:
        "logs/fastqc/{sample}.{unit}.{read}.log"
    wrapper:
        "0.51.2/bio/fastqc"

#
# rule extract_txt:
#     input:
#         "results/qc/fastqc/zip-files/{sample}.{unit}.{read}.fq_fastqc.zip"
#     output:
#         "results/qc/fastqc/txt-files/{sample}.{unit}.{read}.fastqc_data.txt"
#     params:
#         outdir="results/qc/fastqc/txt-files",
#         infile="{sample}.{unit}.{read}_fastqc/fastqc_data.txt"
#     shell:
#         "unzip {input} $(unzip -Z1 {input} | grep 'fastqc_data.txt'$) -d {params.outdir} && mv {params.outdir}/$(unzip -Z1 {input} | grep 'fastqc_data.txt'$) {output} && rmdir {params.outdir}/$(unzip -Z1 {input} | grep 'fastqc_data.txt'$ | cut -d '/' -f1)"
#
#
# files = set()
#
# rule generate_sample_list:
#      input:
#          "results/qc/fastqc/txt-files/{sample}.{unit}.{read}.fastqc_data.txt"
#      output:
#          "results/qc/fastqc/txt-files/test.txt"
#      run:
#
#          samples_list = get_samples_list("results/qc/fastqc/txt-files", ".txt")
#          for i in samples_list:
#              print(i)
#              files.add(i)
#          filepath = os.path.join("results/qc/fastqc/txt-files", "test.txt")
#          f=open("test.txt","w+")
#          f.write("test")
#          f.close()
#          # open("results/qc/fastqc/txt-files/all_files.txt", files)
#
# print(files)

rule multiqc:
    input:
        # files = get_samples_list("results/qc/fastqc/txt-files", ".txt"),
        # expand("results/qc/fastqc/txt-files/{file}.txt", file = files)
        # get_samples_list("results/qc/fastqc/txt-files", ".txt")
        # "results/qc/fastqc/txt-files/{files}.fastqc_data.txt"
        # expand("results/qc/fastqc/txt-files/{sample}.fastqc_data.txt", sample = "{{sample}}")
        # expand("results/qc/fastqc/txt-files/{files}.txt", files=get_samples_list)
        # get_samples_list(directory_path="results/qc/fastqc/txt-files", format="txt")
        # "results/qc/fastqc/txt-files/{files, "."}.fastqc_data.txt"
        # "results/qc/fastqc/txt-files/{sample}.{unit}.{read}.fastqc_data.txt"
        directory("results/qc/fastqc/zip-files")
        # expand("results/qc/fastqc/txt-files/{{sample}}.{{unit}}.{read}.fastqc_data.txt", read=["0", "1", "2"], allow_missing=True)
        # expand("results/qc/fastqc/txt-files/{sample}.{unit}.{read}.fastqc_data.txt", read=["0", "1", "2"], allow_missing=True)
        # dynamic(expand("results/qc/fastqc/txt-files/{{sample}}.{{unit}}.{read}.fastqc_data.txt", read=["0", "1", "2"], allow_missing=True))
    output:
         "results/qc/multiqc/multiqc.html"
    # wildcard_constraints:
    #     sample=units.index.sample,
    #     unit=units.index.unit
    #     files = get_samples_list(directory_path="results/qc/fastqc/txt-files", format="txt")
    #     # file="\w+.\w+"
    params:
        # files = get_samples_list("results/qc/fastqc/txt-files", ".txt")
        ""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "0.51.3/bio/multiqc"
