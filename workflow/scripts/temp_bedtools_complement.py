from snakemake.shell import shell

extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "(bedtools complement"
    " {extra}"
    " -i {snakemake.input.in_file}"
    " -g {snakemake.input.genome}"
    " > {snakemake.output[0]})"
    " {log}"
)
