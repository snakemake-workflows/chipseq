from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell("(preseq lc_extrap {snakemake.params} {snakemake.input[0]} -output {snakemake.output[0]}) {log}")
