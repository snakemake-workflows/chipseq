log = snakemake.log_fmt_shell(stdout=False, stderr=True)

merged_igv = snakemake.input[0]
